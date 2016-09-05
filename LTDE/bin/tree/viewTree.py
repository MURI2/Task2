from __future__ import division
from ete3 import Tree
import os

mydir = os.path.expanduser("~/github/Task2/LTDE")

t = Tree('/Users/WRShoemaker/github/Task2/LTDE/data/Tree/RAxML_MajorityRuleExtendedConsensusTree_test.T18', format=0)

print t
