import os,sys
from optparse import OptionParser

def setupparser(myparser) :
  myparser.add_option("-d", action="store_true", dest="debug", default = False)
  myparser.add_option("-v", action="store_true", dest="verbose", default = False)
  return myparser
