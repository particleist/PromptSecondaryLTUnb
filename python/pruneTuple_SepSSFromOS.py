import os,sys
import ROOT
from ROOT import *
from math import *
from array import *

arraysforouttree = []

def puttrackinouttree(track_P) :
  arraysforouttree[1][0] = track_P.Angle(D0_P)
  arraysforouttree[2][0] = track_P.Angle(K_P)
  arraysforouttree[3][0] = track_P.Angle(pi_P)

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
parser.add_option("--infile", action="store", dest="inputfile", default = "../data/mcd02kpi_tracks5_merged.root")
parser.add_option("--intreename", action="store", dest="intreename", default = "DecayTree")
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

# Get the input tree
infile = TFile(options.inputfile)
intree = infile.Get(options.intreename)
numentries_intree = intree.GetEntries()

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
for i in range(0,3) :
  arraysforouttree.append(array('f',[0]))
outtree.Branch('nTrack',arraysforouttree[0],'nTrack/I')
outtree.Branch('track_angletod',arraysforouttree[1],'track_angletod[nTrack]/F')
outtree.Branch('track_angletokaon',arraysforouttree[2],'track_angletokaon[nTrack]/F')
outtree.Branch('track_angletopion',arraysforouttree[3],'track_angletopion[nTrack]/F')

for entry in range(0,numentries_intree) : 
  intree.GetEntry(entry)
  D0_P = TVector3(intree.D0_PX,intree.D0_PY,intree.D0_PZ)
  K_P  = TVector3(intree.K_PX,intree.K_PY,intree.K_PZ)
  pi_P = TVector3(intree.pi_PX,intree.pi_PY,intree.pi_PZ)
  
  array_hpt_px =  intree.__getattr__('hpt_'+varname+'_px')
  array_hpt_py =  intree.__getattr__('hpt_'+varname+'_py')
  array_hpt_pz =  intree.__getattr__('hpt_'+varname+'_pz')

  array_vdchi2_px =  intree.__getattr__('vdchi2_'+varname+'_px')
  array_vdchi2_py =  intree.__getattr__('vdchi2_'+varname+'_py')
  array_vdchi2_pz =  intree.__getattr__('vdchi2_'+varname+'_pz')

  for track in range(0,array_hpt_px.__len__() ) :
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_hpt_px[track],array_hpt_py[track],array_hpt_pz[track])
    # This bit of dirty coding can be removed once we get rid of  
    # the D0 children in this array in the input ntuple
    if track_P.Angle(K_P) < 0.001 :
      continue
    puttrackinouttree(track_P)
    outtree.Fill()
  for track in range(0,array_vdchi2_px.__len__() ) :
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_vdchi2_px[track],array_vdchi2_py[track],array_vdchi2_pz[track])
    # This bit of dirty coding can be removed once we get rid of  
    # the D0 children in this array in the input ntuple
    if track_P.Angle(K_P) < 0.001 :
      continue
    puttrackinouttree(track_P)
    outtree.Fill()

outfile.Write()
outfile.Close()

