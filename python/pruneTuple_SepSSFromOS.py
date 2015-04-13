import os,sys
import ROOT
from ROOT import *
from math import *
from array import *

arraysforouttree = []
arraysforouttree_doca = []

def puttrackinouttree(track_P,D_P,ndchildren) :
  for i in range(0,ndchildren+1) :
    arraysforouttree[i+1][0] = track_P.Angle(D_P[i])

def distance(M0, V0, M1, V1):
    """ Distance between line 0 passing via M0, of direction
    V0 and line 1 passing through M1, of direction V1 """
    dist = -1
    if V0.Cross(V1).Mag() == 0:
        # Same vector, we have the distance between
        # a point and a line
        M = M1 - M0
        dist = M.Cross(V0).Mag() / V0.Mag()
    else:
        # general case, using teh produit mixte
        M = M1 - M0
        x  = abs(M.Cross(V0) * V1)
        y = (V0.Cross(V1)).Mag()
        dist = x / y
    return dist

def puttrackdocainouttree(track_P, track_FS, D_P, D_FS, ndchildren) :
  for i in range(0,ndchildren+1) :
    arraysforouttree_doca[i+1][0] = distance(track_FS, track_P, D_FS[i], D_P[i])


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
print "Processing number of events: %d" % numentries_intree

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

arraysforouttree_doca.append(array('i',[0]))
for i in range(0,ndchildren+1) :
  arraysforouttree_doca.append(array('f',[0]))


outtree.Branch('nTrack',arraysforouttree[0],'nTrack/I')
outtree.Branch('nTrack',arraysforouttree_doca[0],'nTrack/I')
outtree.Branch('track_angletod',arraysforouttree[1],'track_angletod[nTrack]/F')
outtree.Branch('track_docatod',arraysforouttree_doca[1],'track_docatod[nTrack]/F')
for i in range(1,ndchildren+1) :
  outtree.Branch('track_angletochild'+str(i),arraysforouttree[i+1],'track_angletochild'+str(i)+'[nTrack]/F')
  outtree.Branch('track_docatochild'+str(i),arraysforouttree_doca[i+1],'track_docatochild'+str(i)+'[nTrack]/F')


for entry in range(0,numentries_intree) : 
  # Iterating on tuple entries
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

  # Now preparing the data with the first state position
  D_FS = []
  D_FS.append(TVector3(intree.D_FSX,intree.D_FSY,intree.D_FSZ))
  for i in range(1,ndchildren+1) :
    D_FS.append(TVector3(intree.__getattr__('DChild'+str(i)+'_FSX'),
                        intree.__getattr__('DChild'+str(i)+'_FSY'),
                        intree.__getattr__('DChild'+str(i)+'_FSZ')))

  array_hpt_fsx =  intree.__getattr__('hpt_'+varname+'_FSX')
  array_hpt_fsy =  intree.__getattr__('hpt_'+varname+'_FSY')
  array_hpt_fsz =  intree.__getattr__('hpt_'+varname+'_FSZ')

  array_vdchi2_fsx =  intree.__getattr__('vdchi2_'+varname+'_FSX')
  array_vdchi2_fsy =  intree.__getattr__('vdchi2_'+varname+'_FSY')
  array_vdchi2_fsz =  intree.__getattr__('vdchi2_'+varname+'_FSZ')

  for track in range(0,array_hpt_px.__len__() ) :
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_hpt_px[track],array_hpt_py[track],array_hpt_pz[track])
    track_FS = TVector3(array_hpt_px[track],array_hpt_py[track],array_hpt_pz[track])
    #print track_P, "=", track_FS, "=",  D_P, "=", D_FS, "=", ndchildren
    puttrackdocainouttree(track_P, track_FS, D_P, D_FS, ndchildren)
    outtree.Fill()
  for track in range(0,array_vdchi2_px.__len__() ) :
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_vdchi2_px[track],array_vdchi2_py[track],array_vdchi2_pz[track])
    track_FS = TVector3(array_vdchi2_px[track],array_vdchi2_py[track],array_vdchi2_pz[track])
    puttrackdocainouttree(track_P, track_FS, D_P, D_FS, ndchildren)
    outtree.Fill()


outfile.Write()
outfile.Close()

