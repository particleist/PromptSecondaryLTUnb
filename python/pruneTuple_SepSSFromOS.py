import os,sys
import ROOT
from ROOT import *
from math import *
from array import *

arraysforouttree = []

def puttrackinouttree(track_P,D_P,ndchildren) :
  for i in range(0,ndchildren+1) :
    arraysforouttree[i+1][0] = track_P.Angle(D_P[i])

'''
This script prunes the larger tuple in the "data" directory into two smaller "signal"
and "background" tuples which can be used to train the classifiers, whether
TMVA or scikit-learn or whatever, to separate tracks associated with the 
same side of the decay, and tracks associated with the other side of the decay.
'''

from optparse import OptionParser
parser = OptionParser()
from commonoptions import setupparser
parser = setupparser(parser)
parser.add_option("--infile", action="store", dest="inputfile", default = "../data/mcd02kpi_tracks6_merged.root")
parser.add_option("--intreename", action="store", dest="intreename", default = "DecayTree")
parser.add_option("--trainfraction", action="store", dest="trainfraction", default = "0.1")
parser.add_option("--ndchildren", action="store", dest="ndchildren", default = "2")
parser.add_option("--writeos", action="store_true", dest="writeos", default = False)
parser.add_option("--writess", action="store_true", dest="writess", default = False)
parser.add_option("--outfileos", action="store", dest="outfileos", default = "../data/mcd02kpi_os_forsepssfromos.root")
parser.add_option("--outfiless", action="store", dest="outfiless", default = "../data/mcd02kpi_ss_forsepssfromos.root")
(options, args) = parser.parse_args()

# First of all check, thanks to ROOT's file management we can write the opposite or the same side but not both
if not options.writeos and not options.writess :
  print "You have to either process the opposite side or the same side, exiting."
  sys.exit(1)
if options.writeos and options.writess :
  print "You cannot process both the opposite and same side at the same time, exiting."
  sys.exit(1)
# Set up some common variables based on what we are doing
outfilename = options.outfileos if options.writeos else options.outfiless
varname     = "other" if options.writeos else "same"
try : 
  ndchildren = int(options.ndchildren)
except : 
  print "The number of D children must be an integer, please try again."
  sys.exit(1)

# Get the input tree
infile = TFile(options.inputfile)
intree = infile.Get(options.intreename)
numentries_intree = int(intree.GetEntries()*float(options.trainfraction))

'''
Now we prepare the output trees.

For now we are going to put just "physics inspired" variables in there : the angle
between the track and the D, as well as the angle between the track and each of the D
children. For this we need to know the number of D children and we need to have a
naming convention eventually.
'''

outfile = TFile(outfilename,"recreate")
outtree = TTree(options.intreename,options.intreename)
arraysforouttree.append(array('i',[0]))
for i in range(0,ndchildren+1) :
  arraysforouttree.append(array('f',[0]))
outtree.Branch('nTrack',arraysforouttree[0],'nTrack/I')
outtree.Branch('track_angletod',arraysforouttree[1],'track_angletod[nTrack]/F')
for i in range(1,ndchildren+1) :
  outtree.Branch('track_angletochild'+str(i),arraysforouttree[i+1],'track_angletochild'+str(i)+'[nTrack]/F')

for entry in range(0,numentries_intree) : 
  intree.GetEntry(entry)
  D_P = []
  D_P.append(TVector3(intree.D_PX,intree.D_PY,intree.D_PZ))
  for i in range(1,ndchildren+1) :
    D_P.append(TVector3(intree.__getattr__('DChild'+str(i)+'_PX'),
                        intree.__getattr__('DChild'+str(i)+'_PY'),
                        intree.__getattr__('DChild'+str(i)+'_PZ')))
 
  array_hpt_px =  intree.__getattr__('hpt_'+varname+'_px')
  array_hpt_py =  intree.__getattr__('hpt_'+varname+'_py')
  array_hpt_pz =  intree.__getattr__('hpt_'+varname+'_pz')

  array_vdchi2_px =  intree.__getattr__('vdchi2_'+varname+'_px')
  array_vdchi2_py =  intree.__getattr__('vdchi2_'+varname+'_py')
  array_vdchi2_pz =  intree.__getattr__('vdchi2_'+varname+'_pz')

  for track in range(0,array_hpt_px.__len__() ) :
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_hpt_px[track],array_hpt_py[track],array_hpt_pz[track])
    puttrackinouttree(track_P,D_P,ndchildren)
    outtree.Fill()
  for track in range(0,array_vdchi2_px.__len__() ) :
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_vdchi2_px[track],array_vdchi2_py[track],array_vdchi2_pz[track])
    puttrackinouttree(track_P,D_P,ndchildren)
    outtree.Fill()

outfile.Write()
outfile.Close()

