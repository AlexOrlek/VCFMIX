#!/usr/bin/env python
import os
from pathlib import Path
import pytest
import urllib.request

from src.vcfScan import BinomialTest, vcfScan, lineageScan

SOURCE_DIR = os.path.abspath(__file__)
test_vcf_file = os.path.abspath(os.path.normpath(os.path.join(SOURCE_DIR, '..', '..', 'data', 'testdata', '52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')))

# first download test data (vcf file)
if not os.path.exists(test_vcf_file):
    url = 'https://ora.ox.ac.uk/objects/uuid:5e4ec1f8-e212-47db-8910-161a303a0757/download_file?file_format=x-tar&safe_filename=52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz&type_of_work=Dataset'
    urllib.request.urlretrieve(url, test_vcf_file)


def test_BinomialTest():
    """ test binomial test computation """
    bt = BinomialTest(0.001)

    retVal = bt.compute(0, 0)
    assert retVal == (None, None)

    retVal = bt.compute(1, 1)
    assert retVal == (1, 0)

    retVal = bt.compute(0, 1)
    assert retVal == (1, 0)


def test_vcfScan1():
    """ tests definition of regions """
    v = vcfScan()
    v.add_roi('One', set([1, 2, 3]))
    assert v.roi2psn == {'One': set([1, 2, 3])}
    assert v.psn2roi == {1: {'One'}, 2: {'One'}, 3: {'One'}}

    v.add_roi('Two', set([2, 3, 4]))
    assert v.roi2psn == {'One': set([1, 2, 3]), 'Two': set([2, 3, 4])}
    assert v.psn2roi == {1: set(['One']), 2: set(['One', 'Two']), 3: set(['One', 'Two']), 4: set(['Two'])}

    with pytest.raises(ValueError):
        v.add_roi('Not allowed', set([0]))


def test_vcfScan2():
    """ tests reading from a region when none is specified """
    v = vcfScan()
    v.add_roi('One', set([]))

    inputfile = test_vcf_file
    assert os.path.exists(inputfile), 'Input file does not exist.  Please see README.  You may need to install test data.'

    v.parse(vcffile=inputfile)
    assert len(v.bases.index) == 0


def test_vcfScan_3():
    """ tests reading when the info tag does not exist """
    v = vcfScan(infotag='missing')
    inputfile = test_vcf_file
    assert os.path.exists(inputfile), 'Input file does not exist.  Please see README.  You may need to install test data.'
    v.add_roi('One', [1, 2, 3])

    with pytest.raises(KeyError):
        v.parse(vcffile=inputfile)


def test_vcfScan_4():
    """ tests reading from a region """
    v = vcfScan()
    v.add_roi('One', set([1, 2, 3]))
    v.add_roi('Two', set([2, 3, 4]))

    inputfile = test_vcf_file
    assert os.path.exists(inputfile), 'Input file does not exist.  Please see README.  You may need to install test data.'
    v.parse(vcffile=inputfile)
    assert len(v.bases.index) == 6


def test_lineageScan():
    """ tests Loading branch information for deep branches """
    v = lineageScan()
    inputfile = test_vcf_file
    assert os.path.exists(inputfile), 'Input file does not exist.  Please see README.  You may need to install test data.'

    res = v.parse(vcffile=inputfile, guid='528')
    assert res is None

    assert len(v.region_stats.index) == 64

    # check file export works
    targetdir = os.path.join('..', 'misc', 'unitTest_tmp')
    Path(targetdir).mkdir(parents=True, exist_ok=True)
    outputfile = os.path.join(targetdir, '528.txt')
    if os.path.exists(outputfile):
        os.unlink(outputfile)

    v.region_stats.to_csv(outputfile)
    assert os.path.exists(outputfile), 'outputfile does not exist.'

    # compute summary with stored data
    summary1 = v.f_statistics()

    assert isinstance(summary1, dict)
    assert set(summary1.keys()) == set(['mixture_quality', 'F2', 'F47'])

    # compute summary with persisted csv data
    summary2 = v.f_statistics(outputfile)
    assert summary1 == summary2
