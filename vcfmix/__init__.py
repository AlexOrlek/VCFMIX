import importlib.metadata

__version__ = importlib.metadata.version("vcfmix")

from .vcfScan import FastaMixtureMarker, BinomialTest, vcfScan, lineageScan
