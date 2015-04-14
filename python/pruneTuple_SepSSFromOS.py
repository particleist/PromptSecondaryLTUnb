import os,sys
import ROOT
from ROOT import *
from math import *
from array import *

arraysforouttree = []
arraysforouttree_doca = []
arraysforouttree_ddist = []
arraysforttype = []

def puttrackinouttree(track_P,D_P,ndchildren) :
  for i in range(0,ndchildren+1) :
    arraysforouttree[i+1][0] = track_P.Angle(D_P[i])

def distance(M0, V0, M1, V1):
    """ Distance between line 0 passing via M0, of direction
    V0 and line 1 passing through M1, of direction V1 """
    dist = -1
    if V0.Cross(V1).Mag() == 0:
        # Colinear vectors, we have the distance between
        # a point and a line
        M = M1 - M0
        dist = M.Cross(V0).Mag() / V0.Mag()
    else:
        # general case
        M = M1 - M0
        x  = abs(M.Cross(V0) * V1)
        y = (V0.Cross(V1)).Mag()
        dist = x / y
    return dist

def distanceToPoint(M0, V0, M1):
    """ Distance between line 0 passing via M0, of direction
    V0 and point M1 """
    M = M1 - M0
    dist = M.Cross(V0).Mag() / V0.Mag()
    return dist

def puttrackdocainouttree(track_P, track_FS, D_P, D_FS, ndchildren) :
  for i in range(0,ndchildren+1) :
    arraysforouttree_doca[i+1][0] = distance(track_FS, track_P, D_FS[i], D_P[i])

def puttrackddistinouttree(track_P, track_FS, Vertex) :
  arraysforouttree_ddist[0][0] = distanceToPoint(track_FS, track_P, Vertex)

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


# Preparing branches for the D and its children
#########################################################################
arraysforouttree.append(array('i',[0]))
for i in range(0,ndchildren+1) :
  arraysforouttree.append(array('f',[0]))

arraysforouttree_doca.append(array('i',[0]))
for i in range(0,ndchildren+1) :
  arraysforouttree_doca.append(array('f',[0]))

arraysforouttree_ddist.append(array('f',[0]))
arraysforttype.append(array('i',[0]))

outtree.Branch('tracktype', arraysforttype[0],'tracktype/I')
outtree.Branch('nTrack',arraysforouttree[0],'nTrack/I')
outtree.Branch('track_angletod',arraysforouttree[1],'track_angletod[nTrack]/F')
outtree.Branch('track_docatod',arraysforouttree_doca[1],'track_docatod[nTrack]/F')
outtree.Branch('track_devdist',arraysforouttree_ddist[0],'track_devdist[nTrack]/F')

for i in range(1,ndchildren+1) :
  outtree.Branch('track_angletochild'+str(i),arraysforouttree[i+1],'track_angletochild'+str(i)+'[nTrack]/F')
  outtree.Branch('track_docatochild'+str(i),arraysforouttree_doca[i+1],'track_docatochild'+str(i)+'[nTrack]/F')

# Iterating on tuple entries
#########################################################################
for entry in range(0,numentries_intree) : 
  intree.GetEntry(entry)

  # Data related of the momentum of the D
  #########################################################################
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

  array_hpt_fsx =  intree.__getattr__('hpt_'+varname+'_FSX')
  array_hpt_fsy =  intree.__getattr__('hpt_'+varname+'_FSY')
  array_hpt_fsz =  intree.__getattr__('hpt_'+varname+'_FSZ')

  array_vdchi2_fsx =  intree.__getattr__('vdchi2_'+varname+'_FSX')
  array_vdchi2_fsy =  intree.__getattr__('vdchi2_'+varname+'_FSY')
  array_vdchi2_fsz =  intree.__getattr__('vdchi2_'+varname+'_FSZ')

  D_FS = []
  D_FS.append(TVector3(intree.D_FSX,intree.D_FSY,intree.D_FSZ))
  for i in range(1,ndchildren+1) :
    D_FS.append(TVector3(intree.__getattr__('DChild'+str(i)+'_FSX'),
                        intree.__getattr__('DChild'+str(i)+'_FSY'),
                        intree.__getattr__('DChild'+str(i)+'_FSZ')))

  D_EV = TVector3(intree.D_VX,intree.D_VY,intree.D_VZ)

  for track in range(0,array_hpt_px.__len__() ) :
    arraysforttype[0][0] = 1
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_hpt_px[track],array_hpt_py[track],array_hpt_pz[track])
    puttrackinouttree(track_P,D_P,ndchildren)

    track_FS = TVector3(array_hpt_px[track],array_hpt_py[track],array_hpt_pz[track])
    puttrackdocainouttree(track_P, track_FS, D_P, D_FS, ndchildren)

    puttrackddistinouttree(track_P, track_FS, D_EV)

    outtree.Fill()

  for track in range(0,array_vdchi2_px.__len__() ) :
    arraysforttype[0][0] = 2
    arraysforouttree[0][0] = 1
    track_P = TVector3(array_vdchi2_px[track],array_vdchi2_py[track],array_vdchi2_pz[track])
    puttrackinouttree(track_P,D_P,ndchildren)

    arraysforouttree_doca[0][0] = 1
    track_P = TVector3(array_vdchi2_px[track],array_vdchi2_py[track],array_vdchi2_pz[track])
    puttrackdocainouttree(track_P, track_FS, D_P, D_FS, ndchildren)

    puttrackddistinouttree(track_P, track_FS, D_EV)

    outtree.Fill()


outfile.Write()
outfile.Close()

