#!/usr/bin/env python3
"""
classes to parse VCF files extracting high-quality bases from a pre-defined set of positions.

Applications include detection of mixtures of TB genomes at lineage defining sites, or sites varying within a genome.
These classes can create a persistent, indexed extract of a subset of data.

For the regionScan_from_genbank class, which computes variation over
regions of a genbank file, please see https://github.com/davidhwyllie/adaptivemasking.


Note that in general, some methods within these classes require a sequence identifier to be provided.
These are referred to as 'guids' and we hav always used globally unique identifiers (guids) as sequence identifiers.
However, the requirement to actually use a guid is not absolute, and is not enforced.
Other character strings identifying a sequence should also work, subject to:
* not being more than 36 characters
* being valid as part of a file name
"""

import os
import gzip
from collections import deque
from scipy import stats
import numpy as np
import math
import logging
import warnings
import pandas as pd
import tables

from Bio import SeqIO

SOURCE_DIR = os.path.abspath(__file__)


class FastaMixtureMarker():
    """ writes codes reflecting mixed base calls into a fasta file """

    def __init__(self, expectedErrorRate, mlp_cutoff, clustering_cutoff=None, min_maf=0):
        """ writes mixed base calls into fasta files.
            expectedErrorRate: the minor variant frequency expected by chance, given mapper performance.  Typically 0.001 if Q30 cutoffs used when mapping.
            mlp_cutoff : the -log P cutoff used for mixed base selection.
            min_maf: the minimum minor variant frequency reported
            clustering_cutoff: call bases N not mixed if within clustering_cutoff of another mixed base
                rationale is that with some sequencer/mapper/basecalling combinations, such clustered mixed bases may reflect an underlying genomic feature, such as an indel,
                rather than a genuine mix of sequences.
        """
        self.expectedErrorRate = expectedErrorRate
        self.mlp_cutoff = mlp_cutoff
        self.min_maf = min_maf
        self.clustering_cutoff = clustering_cutoff
        self.bt = BinomialTest(expectedErrorRate)
        # iupac codes from https://www.bioinformatics.org/sms/iupac.html
        # capitalisation indicates which way round the base frequencies are.
        # if lower case, the first base alphabetically is more common;
        # if upper case, the second is more common.
        self.iupac = {
            'AG': 'r', 'GA': 'R',
            'AT': 'w', 'TA': 'W',
            'CT': 'y', 'TC': 'Y',
            'AC': 'm', 'CA': 'M',
            'CG': 's', 'GC': 'S',
            'GT': 'k', 'TG': 'K'
        }

    def mark_mixed(self, seq_file, mixed_bases_file):
        """ annotations fasta with mixed bases, using IUPAC codes.

        Inputs:
            seq_file: a fasta file.
            mixed_bases_file: a csv file containing mixed bases, e.g. as generated by v.parse() followed by v.bases.to_csv(filename), where v is a vcfScan object.

        Outputs:
            a tuple, consisting of a Pandas dataframe containing the positions marked, and a string comprising the fasta string with IUPAC codes added, OR
            a tuple consisting of (None,e) where e is an exception .
        """

        # read fasta
        if seq_file.endswith('.gz'):
            try:
                with gzip.open(seq_file, 'rt') as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        seq = list(record.seq)
            except Exception as e:
                return None, e
        else:
            try:
                with open(seq_file, 'r') as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        seq = list(record.seq)
            except Exception as e:
                return None, e
        # read outputfile, if there are any columns
        try:
            df = pd.read_csv(mixed_bases_file, index_col='pos')
        except pd.errors.EmptyDataError:
            return None, ''.join(seq)

        variants_to_update = {}
        for ix in df.index:

            maf = df.loc[ix, 'maf']

            if maf >= self.min_maf:
                mlp = df.loc[ix, 'mlp']

                if np.isnan(mlp):
                    # compute it
                    p_value, mlp = self.bt.compute(df.loc[ix, 'nonmajor'], df.loc[ix, 'depth'])

                if mlp >= self.mlp_cutoff:		# the cutoff
                    base = pd.DataFrame({'base': ['A', 'C', 'G', 'T'],
                                         'depth': [df.loc[ix, 'base_a'],
                                                   df.loc[ix, 'base_c'],
                                                   df.loc[ix, 'base_g'],
                                                   df.loc[ix, 'base_t']
                                                   ]
                                         })

                    base = base.sort_values(by='depth', ascending=False)
                    top2 = ''.join(base.head(2)['base'].tolist())
                    variants_to_update[ix - 1] = {'pos': ix - 1, 'base': self.iupac[top2]}

        # eliminate clustered Ms if self.clustering_cutoff is not None:
        bases = sorted(variants_to_update.keys())
        if self.clustering_cutoff is not None:
            for i in range(1, len(bases) - 1):
                if bases[i] - bases[i - 1] <= self.clustering_cutoff:
                    variants_to_update[bases[i]]['base'] = 'N'
                if bases[i + 1] - bases[i] <= self.clustering_cutoff:
                    variants_to_update[bases[i]]['base'] = 'N'

        for i in range(1, len(bases) - 1):
            seq[bases[i]] = variants_to_update[bases[i]]['base']

        df = pd.DataFrame.from_dict(variants_to_update, orient='index')
        if len(df.index) > 0 and len(df.columns.values.tolist()) > 0:
            if 'base' in df.columns.values.tolist() and 'pos' in df.columns.values.tolist():
                df = df.query("not base=='N'")
                df = df.query("pos>0")		# don't count the first base
        return df, ''.join(seq)


class BinomialTest():
    """ wrapper round stats.binom_test, testing the significance of a particular minor variant count
    given a particular depth, and expected error rate.

    Stores the results of binomial tests in a dictionary, which results in faster computations as the
    computation is only done once, then stored. """

    def __init__(self, expectedErrorRate):
        self.expectedErrorRate = expectedErrorRate
        self.p_values = {}
        if not isinstance(expectedErrorRate, float):
            raise TypeError("expectedErrorRate supplied is {0}; this must but a float, not a {1}".format(expectedErrorRate, type(expectedErrorRate)))

    def compute(self, minor_variant_count, depth):
        """ compute binomial test.  returns p value and minus log p. """
        if depth == 0:
            return None, None
        if minor_variant_count == depth:
            return 1, 0

        key = "{0},{1}".format(minor_variant_count, depth)
        try:
            p_value = self.p_values[key]
        except KeyError:
            # not precomputed
            self.p_values[key] = stats.binom_test(x=minor_variant_count, n=depth, p=self.expectedErrorRate)   # do the test if any variation

        p_value = self.p_values[key]
        if p_value == 0:
            mlp = 250        # code minus log p as 250 if p value is recorded as 0 in float format
        elif p_value is not None:
            mlp = -math.log(p_value, 10)
        elif p_value is None:
            mlp = None
        return p_value, mlp


class vcfScan():
    """ parses a VCF file, extracting high quality calls at pre-defined positions.

            Note that at present this doesn't support BCF files.
            To enable this, we'd need to modify the ._parse function
            to use pysam https://pypi.python.org/pypi/pysam.

    """

    def __init__(self,
                 expectedErrorRate=0.001,
                 infotag='BaseCounts4',
                 report_minimum_maf=0,
                 compute_pvalue=True):
        """ creates a vcfScan object.

        Args:
            expectedErrorRate: a floating point number indicating the expected per-base error rate.
            infotag: a string indicating the vcf/bcf INFO tag from which to extract the four high quality depths corresponding to A,C,G,T.
                If infotag == 'AD', it is assumed it is in the format exported by *samtools mpileup*.
                If infotag == 'auto', it will look for AD tags, then BaseCounts4 tags, and if both are missing, fail.
                Otherwise, it assumes it is in the format generated by GATK VariantAnnotator.
            report_minimum_maf: a float between 0 and 1.  Bases with mixed allele frequency (maf)
                less than report_minimum_maf will not be reported.  If report_minimum_maf > 0,
                bases with zero depth will not be reported.  If low-frequency minor variants (e.g. <5%)
                are not of interest, setting report_minimum_maf markedly increases speed and reduces memory requirements.
            compute_pvalue: if True, performs an exact binomial test per base, comparing the observed minor variant
                frequency with expectedErrorRate

        Returns:
            None

        """
        self.roi2psn = dict()       # region of interest -> genomic position
        self.psn2roi = dict()       # genomic position -> region of interest
        self.infotag = infotag
        self.fieldtag = None
        self.report_minimum_maf = report_minimum_maf
        self.compute_pvalue = compute_pvalue
        self.expectedErrorRate = expectedErrorRate
        self.bt = BinomialTest(self.expectedErrorRate)

    def add_roi(self, roi_name, roi_positions):
        """ adds a region of interest, which is a set of genomic positions for which the
        variation should be extracted.

        self.roi2psn and self.psn2roi allow in memory lookup between positions and regions, and v/v;
        add_roi creates entries in roi2psn and psn2roi for all roi_names and roi_positions.

        Note that the roi_positions must be 1-indexed.  A value error is raised if a zero position is added.

        Args:
            roi_name: the name of the region of interest, example 'gene3'
            roi_position: a list containing the one indexed positions  of the bases in roi_name

        Returns:
            None
        """

        try:
            self.roi2psn[roi_name] = set([])
        except KeyError:
            # already exists
            pass
        for roi_position in roi_positions:
            if roi_position == 0:
                raise ValueError("Positions supplied must be 1, not zero indexed")

            self.roi2psn[roi_name].add(roi_position)
            if roi_position in self.psn2roi.keys():
                self.psn2roi[roi_position].add(roi_name)
            else:
                self.psn2roi[roi_position] = set([roi_name])

    def persist(self, outputfile, mode='w'):
        """ persists any bases examined (i.e. which are part of rois) to an indexed hdf5 file.

        For M. tuberculosis, which has a 4.4e6 nt genome, this takes up about 52MB if
        all bases are part of ROIs (as, for example, implemented by regionScan).

        The bases examined are indexed by ROI and position,
        allowing near instantaneous access from on-disc stores.

        The HDF store access is implemented via Pandas and PyTables.
        Each matrix stored is associated with a key.
        The key used is self.guid, which is set by self.parse().
        If guid is not set, an error is raised.

        By default, any existing HDF file will be overwritten (mode 'w').
        To append data to an existing hdf file, use mode 'a'.

        Args:
            outputfile: the outputfile name
            mode: 'a' to append to an existing file; 'w' to overwrite.

        Returns:
            None

        """

        self._persist(self.bases, outputfile, mode)

    def _persist(self, df, outputfile, mode='w'):
        """ persists df to an indexed hdf5 file.

        Args:
            df: the data frame to export
            outputfile: the outputfile name
            mode: 'a' to append to an existing file; 'w' to overwrite.

        Returns:
            None

        """

        if self.guid is None:
            raise ValueError("Cannot write hdf file with a null table key.  You must set the guid when you .parse() the vcf file")
        if df is None:
            raise ValueError("No data to write; None passed.")

        warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)		# guids are not valid python names, but this doesn't matter
        df.to_hdf(outputfile,
                  key=self.guid,
                  mode=mode,
                  format='t',
                  complib='blosc:blosclz',
                  data_columns=['roi_name', 'pos'],
                  complevel=9)

    def parse(self, vcffile, guid=None):
        """ parses a vcffile.
        stores a pandas data frame, with one row per roipsn/roiname combination, in self.bases.
        You must provide a guid is you wish to persist the object using the .persist method.
        Wrapper around ._parse().

        Args:
            vcffile: the vcf file to read
            guid: a guid identifier for the parsed object; required only if using .persist() to store the parsed object.

        Returns:
            True, if the parse succeeded
            False, if it did not, for example due to file truncation
        """
        self.guid = guid
        self._parse(vcffile)

    def _parse(self, vcffile):
        """ parses a vcffile.
        stores a pandas data frame, with one row per roipsn/roiname combination, in self.bases

        Args:
            vcffile: the vcf file to parse

        Returns:
            None

        Notes:
            Output is stored in self.bases
        """

        # set up variable for storing output
        resDict = {}
        nAdded = 0
        self.region_stats = None
        self.bases = None
        warning_emitted = False

        # transparently handle vcf.gz files.
        if vcffile.endswith('.gz'):
            f = gzip.open(vcffile, "rb")
        else:
            f = open(vcffile, "rt")

        # precompute a sorted list of positions to look for
        # in an efficient data structure
        sought_psns = deque(sorted(self.psn2roi.keys()))

        try:
            sought_now = sought_psns.popleft()
            if sought_now == 0:
                warnings.warn("Asked to estimate mixtures for base 0.  Positions should be 1 -indexed")
            # iterate over the vcf file
            for line in f:

                if not isinstance(line, str):
                    line = line.decode()		# turn it into a string

                if line[0] == "#":
                    continue  # it is a comment; go to next line;

                if "INDEL" in line:
                    continue  # this is not needed because we're not dealing with indels here; next line;

                # parse the line.
                chrom, pos, varID, ref, alts, score, filterx, infos, fields, sampleInfo = line.strip().split()
                pos = int(pos)

                # the current position (pos) should be <= sought_now
                # if all bases in the vcf are called.
                # if they are not, then we need to 'catch up' and find the next
                # sought_now position after or at the current vcf scan position, pos
                if not pos <= sought_now:
                    if warning_emitted is False:
                        logging.warn("Note: not all positions are called in vcf file: gap observed nr bases {1}..{0}; adjusting scan.\
                            Results should not be affected. Subsequent similar warnings will not be shown.".format(pos, sought_now))
                        warning_emitted = True
                    while sought_now <= pos:
                        try:
                            sought_now = sought_psns.popleft()
                        except IndexError:		# no more positions defined for selection; this is allowed
                            # we're out of positions
                            sought_now = 1e12 	# bigger than any genome; will never get there

                if pos == sought_now:
                    # we are looking for a position in the vcf file, and we have found it;
                    alts = alts.split(",")
                    alts = [i for i in alts if i in ['A', 'T', 'C', 'G']]
                    infos = dict(item.split("=") for item in infos.split(";"))
                    fields = fields.split(':')
                    sampleInfo = sampleInfo.split(':')

                    # autodetect: if self.infotag is 'auto' and AD tag present, use that
                    if self.infotag == 'auto':
                        if 'AD' in fields:
                            self.fieldtag = 'AD'
                            self.infotag = None
                        elif 'AD' in infos.keys():
                            self.infotag = 'AD'
                        elif 'BaseCounts4' in infos.keys():  # BaseCounts4 is the high-quality per-base depth from GATK VariantAnnotator
                            self.infotag = 'BaseCounts4'
                        else:
                            raise KeyError("auto detection of infotag required, but neither AD nor BaseCounts4 tags found ")

                    # parse self.infotag/self.fieldtag to extract baseCounts
                    if self.infotag is not None:  # infotag is either AD or BaseCounts4
                        try:
                            baseCounts = list(map(int, infos[self.infotag].split(",")))  # get frequencies of high quality bases
                        except KeyError:
                            raise KeyError("Expected a tag {0} in the 'info' component of the call file, but it was not there.  Keys present are: {1}".format(self.infotag, infos.keys()))
                        except ValueError:	 # very rare - likely reflects vcf corruption
                            warnings.warn("Integer conversion failed at VCF row {0}: was applied to {1}; line is: {2}.".format(pos, infos[self.infotag], line))
                            baseCounts = [0, 0, 0, 0]  # assign zero basecounts

                    if self.fieldtag == 'AD':
                        try:
                            baseCounts = list(map(int, sampleInfo[fields.index(self.fieldtag)].split(",")))
                        except KeyError:
                            raise KeyError("Expected a tag {0} in the 'fields' component of the call file, but it was not there.  Keys present are: {1}".format(self.fieldtag, fields))

                    if self.infotag == 'AD' or self.fieldtag == 'AD':  # then this is a vcf file made by samtools mpileup with -t INFO/AD or -t AD flag
                        # and we need to extract the data accordingly from the ref and alt columns.
                        # bases which are not mentioned in the ALT are zero
                        basedict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                        # gather basecounts from the AD tag; first element is REF depth, subsequent elements are ALT depths
                        basedict[ref] = baseCounts[0]
                        counter = 0
                        for alt in alts:
                            counter += 1
                            basedict[alt] = baseCounts[counter]
                        # and generate a baseCounts corresponding to that produced by GATK VariantAnnotator
                        baseCounts = [basedict['A'], basedict['C'], basedict['G'], basedict['T']]

                    # extract the baseCounts, and do QC
                    baseFreqs = baseCounts.copy()
                    baseFreqs.sort(reverse=True)
                    if not len(baseFreqs) == 4:
                        raise TypeError("Expected tag {0} to contain 4 depths, but {1} found.  Base = {2}; tag contents are {3}".format(self.infos, len(baseFreqs), pos, baseCounts))
                    depth = sum(baseCounts)

                    # compute probability that the minor variant frequency differs from self.expectedErrorRate from exact binomial test
                    pvalue = None
                    mlp = None
                    if self.compute_pvalue:
                        pvalue, mlp = self.bt.compute(baseFreqs[1] + baseFreqs[2] + baseFreqs[3], depth)

                    # store output in a dictionary
                    if depth > 0:
                        maf = float(baseFreqs[1]) / float(depth)
                    else:
                        maf = None

                    for roi_name in self.psn2roi[sought_now]:
                        report_base = True
                        if maf is None and self.report_minimum_maf > 0:
                            report_base = False
                        if maf is not None:
                            if maf < self.report_minimum_maf:
                                report_base = False

                        if report_base:
                            nAdded += 1
                            resDict[nAdded] = {'roi_name': roi_name, 'pos': pos, 'ref': ref, 'depth': depth,
                                               'base_a': baseCounts[0],
                                               'base_c': baseCounts[1],
                                               'base_g': baseCounts[2],
                                               'base_t': baseCounts[3],
                                               'maf': maf,
                                               'mlp': mlp}

                    # recover the next item to recover
                    try:
                        sought_now = sought_psns.popleft()
                    except IndexError:		# no positions selected
                        break				# all positions have been selected

        except IndexError:		# no positions defined for selection; this is allowed
            pass
        except EOFError:		# file is corrupt
            return False		# failed
        # construct data frame
        self.bases = pd.DataFrame.from_dict(resDict, orient='index')

        # construct summary by region, defined by roi_name
        if len(self.bases.index) > 0:
            r1 = self.bases.groupby(['roi_name'])['depth'].mean().to_frame(name='mean_depth')
            r2 = self.bases.groupby(['roi_name'])['depth'].min().to_frame(name='min_depth')
            r3 = self.bases.groupby(['roi_name'])['depth'].max().to_frame(name='max_depth')

            r4 = self.bases.groupby(['roi_name'])['pos'].min().to_frame(name='start')
            r5 = self.bases.groupby(['roi_name'])['pos'].max().to_frame(name='stop')
            r6 = self.bases.groupby(['roi_name'])['pos'].count().to_frame(name='length')

            # if all mafs are NA, then mean() will fail with a pandas.core.base.DataError
            try:
                r8 = self.bases.groupby(['roi_name'])['maf'].mean().to_frame(name='mean_maf')
            except pd.core.base.DataError:
                r8 = r1.copy()
                r8.columns = ['mean_maf']
                r8['mean_maf'] = None

            # compute total depth
            r9 = self.bases.groupby(['roi_name'])['depth'].sum().to_frame(name='total_depth')

            # compute total_nonmajor_depth
            self.bases['most_common'] = self.bases[['base_a', 'base_c', 'base_g', 'base_t']].max(axis=1)
            self.bases['nonmajor'] = self.bases['depth'] - self.bases['most_common']
            r10 = self.bases.groupby(['roi_name'])['nonmajor'].sum().to_frame(name='total_nonmajor_depth')

            df = pd.concat([r1, r2, r3, r4, r5, r6, r8, r9, r10], axis=1)              # in R,  this is a cbind operation
        else:
            df = None
        self.region_stats = df
        f.close()
        return True


class lineageScan(vcfScan):
    """ parses a vcf file, extracting high-quality bases at lineage defining positions, as
    described by Coll et al. """

    def __init__(self,
                 expectedErrorRate=0.001,
                 lineage_definition_file=os.path.abspath(os.path.normpath(os.path.join(SOURCE_DIR, '..', '..', 'data', 'refdata', 'Coll2014_LinSpeSNPs_final.csv'))),
                 exclusion_position_file=os.path.abspath(os.path.normpath(os.path.join(SOURCE_DIR, '..', '..', 'data', 'refdata', 'exclusion_nt.txt'))),
                 infotag='BaseCounts4'
                 ):
        """ creates a vcfScan object.
        Inputs:
            expectedErrorRate expected error rate in basecalling.  Used in a negative binomial test.
            lineage_definition_file a file containing a classification of lineage defining positions, e.g. Coll et al.
            exclusion_position_file a set of positions not to call.  A default file for H37Rv v2 is provided.

        Outputs:
            self.lineage_defining contains the lineage defining positions read

        """
        self.roi2psn = dict()
        self.psn2roi = dict()
        self.expectedErrorRate = expectedErrorRate
        self.infotag = infotag
        self.fieldtag = None
        self.excluded = set()
        self.report_minimum_maf = 0
        self.compute_pvalue = True
        if exclusion_position_file is not None:
            excluded_positions = pd.read_table(exclusion_position_file, sep=',', header=0)
            self.excluded = set(excluded_positions['pos'].tolist())

        self.lineage_defining = pd.read_table(lineage_definition_file, sep=',', header=0)
        self.lineage_defining['pos'] = self.lineage_defining['position']

        for roi_name in self.lineage_defining['lineage'].unique():
            roi_psns = set(self.lineage_defining[self.lineage_defining['lineage'] == roi_name]['position'].tolist())
            roi_psns = roi_psns - self.excluded		# remove any high variation regions
            self.add_roi(roi_name, roi_psns)
        self.bt = BinomialTest(self.expectedErrorRate)

    def parse(self, vcffile, guid):
        """ parses a vcffile.

        Inputs:
            vcffile: the vcffile to scan

        Outputs:
            self.region_stats contains minor variant calls in the regions defined.
            self.bases are the bases included in the regions of interest (only)

        Returns:
            None

        To export output, after calling parse to
        obj.region_stats.to_csv(filename)

        To generate f2 and f47 statistics, as described in publication, call
        obj.f_statistics()

        """

        self._parse(vcffile)
        self.region_stats['guid'] = guid
        return(None)

    def f_statistics(self, filename=None):
        """ computes F2 and F47 summary statistics.
            If filename is None (default) then it computes them on self.region_stats.
            If the filename is not none, and exists, it assumes that the file is csv data containing region statistics, as
            exported by obj.region_stats.to_csv(filename).  It loads this data, and reports it.
        """

        if filename is not None:  # compute statistics on a stored summary file;
            if os.path.exists(filename):
                self.region_stats = pd.read_csv(filename)
                # sanity check
                existing_columns = set(self.region_stats.columns.values.tolist())
                expected_columns = set(['roi_name', 'mean_depth', 'min_depth', 'max_depth', 'start', 'stop', 'length',
                                        'mean_maf', 'total_depth', 'total_nonmajor_depth', 'guid'])
                if not existing_columns == expected_columns:
                    raise KeyError("Read filename {0} but the data frame had columns {1} not {2} diffs: {3}; {4}".format(
                        filename,
                        existing_columns,
                        expected_columns,
                        existing_columns - expected_columns,
                        expected_columns - existing_columns))
            else:
                raise FileExistsError("Asked to read lineage summary from {0} which does not exist".format(filename))

        if self.region_stats is not None:
            if len(self.region_stats) < 58:		# lineage defining sites were not computed
                return({'mixture_quality': 'bad',
                        'F2': None,
                        'F47': None})

            else:
                sorted_region_stats = self.region_stats.sort_values(by='mean_maf', ascending=False)
                f2_denominator = sum(sorted_region_stats['total_depth'].head(2))
                f47_denominator = sum(sorted_region_stats['total_depth'].tail(47))

                # trap for the situation in which there are no reads, so F2 and F47 can't be computed (divide-by-zero)
                if f2_denominator == 0 or f47_denominator == 0:
                    return({'mixture_quality': 'bad',
                            'F2': None,
                            'F47': None})
                else:
                    f2 = sum(sorted_region_stats['total_nonmajor_depth'].head(2)) / f2_denominator
                    f47 = sum(sorted_region_stats['total_nonmajor_depth'].tail(47)) / f47_denominator

                    return({'mixture_quality': 'OK',
                            'F2': f2,
                            'F47': f47})
