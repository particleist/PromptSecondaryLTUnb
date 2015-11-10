#!/usr/bin/env python
# =============================================================================
import ROOT
from Bender.MainMC import * 
from GaudiKernel.SystemOfUnits import GeV, MeV, mm
import BenderTools.Fill

# =============================================================================
## (optional) logging
# =============================================================================

__author__  = " Ben Couturier ben.couturier@cern.ch "
__date__    = " 2014-02-05 " 
__version__ = " Version $Revision:$ "

from Bender.Logger import getLogger 
if '__main__' == __name__ : logger = getLogger ( 'Bender.Startup' )
else                      : logger = getLogger ( __name__ )
# =============================================================================

## @class TrackFilter
#  Simple class to filter tracks in D0 -> K pi decays
class TrackFilter(AlgoMC):
    """
    Algorithm to filter interesting tracks from the same PV as the foound D0
    """

    ## initialize the algorithm
    def initialize ( self ) :
        """
        Initialize the algorithm
        """
        ## initialize the base
        sc = AlgoMC.initialize ( self ) ## initialize the base
        if sc.isFailure() :
            return sc
        
        ## initialize functionality for enhanced n-tuples 
        sc = self.fill_initialize ()     ## IMPORTANT 

        ## Loading the particle 2 MC associators
        # This is taken from TupleToolMCTruth.cpp
        self._associators = [ self.tool(cpp.IParticle2MCAssociator, tool)
                              for tool in ["DaVinciSmartAssociator",
                                           "MCMatchObjP2MCRelator"]]
        return SUCCESS

    ## Make sure we cleanup in finalize
    def finalize ( self ) :    
        self.fill_finalize  ()                  ## IMPORTANT 
        self._associators = None
        return AlgoMC.finalize ( self )

    ## Util method for looking up associated monte carlo particle
    def _getRelatedMCParticle(self, p):
        """
        Method that looks for the related MC particle using the associators initialized
        in _associators
        """
        mcp = None
        for a in self._associators:
            mcp = a.relatedMCP(p)
            if mcp: break
        return mcp

    ## Find the top of a MC decay for a givent MCParticle
    def _findMCTop(self, mcp):
        """
        Iterate through the parents to find the top decay
        """
        if mcp == None:
            return None
        tmp = mcp
        while tmp.mother() != None:
            tmp = tmp.mother()
        return tmp

    ## Find the top of a MC decay for a givent MCParticle
    def _isFromB(self, mcp):
        """
        Iterate through the parents to find the top decay
        """
        if mcp == None:
            return 0

        isfromB = 0
        tmp = mcp
        while tmp.mother() != None:
            motherId = abs(tmp.mother().particleID().pid())
            if ( motherId >= 500 and motherId < 600 ) or ( motherId >= 5000 and motherId < 6000 ):
                isfromB = 1
                break
            tmp = tmp.mother()
            
        return isfromB

    ## Dumping particle info
    def processParticle(self, p, tuple, prefix, full=False, firstState=False):
        """
        Dump particle info to a tuple
        """
        tuple.column(prefix + "ID", int(ID(p)))
        tuple.column(prefix + "PX", PX(p))
        tuple.column(prefix + "PY", PY(p))
        tuple.column(prefix + "PZ", PZ(p))
        tuple.column(prefix + "PT", PT(p))
        tuple.column(prefix + "E", E(p))
        tuple.column(prefix + "ETA", ETA(p))
        tuple.column(prefix + "PHI", PHI(p))
        tuple.column(prefix + "MM", MM(p) )
        tuple.column(prefix + "VX",  VX(p.endVertex()))
        tuple.column(prefix + "VY",  VY(p.endVertex()))
        tuple.column(prefix + "VZ",  VZ(p.endVertex()))
        tuple.column(prefix + "BPVVD",  BPVVD(p))
        tuple.column(prefix + "BPVVDCHI2",  BPVVDCHI2(p))
        tuple.column(prefix + "BPVIPCHI2", BPVIPCHI2()(p) )
	if firstState:
            fsx = -1e6
            fsy = -1e6
            fsz = -1e6
            if p.proto() != None:
                if p.proto().track() != None:
                    fs = p.proto().track().firstState()
                    if fs != None:
                        fsx = fs.x()
                        fsy = fs.y()
                        fsz = fs.z()
            tuple.column(prefix  + 'FSX', fsx)
            tuple.column(prefix  + 'FSY', fsy)
            tuple.column(prefix  + 'FSZ', fsz)

        if full:
            tuple.column(prefix + "BPVDIRA",  BPVDIRA(p))
            tuple.column(prefix + "BPVLTIME",   BPVLTIME()(p))
            tuple.column(prefix + "BPVLTCHI2",   BPVLTCHI2()(p))
            bpv = self.getRelatedPV(p)
            tuple.column(prefix + "PVX",  VX(bpv))
            tuple.column(prefix + "PVY",  VY(bpv))
            tuple.column(prefix + "PVZ",  VZ(bpv))

    def processParticlePID(self, p, tuple, prefix):
         """
         Dump particle info to a tuple
         """
         tuple.column(prefix + "PIDK", PIDK(p))
         tuple.column(prefix + "PIDp", PIDp(p))
         tuple.column(prefix + "PIDpi", PIDpi(p))
         tuple.column(prefix + "PIDmu", PIDmu(p))
         tuple.column(prefix + "ID", int(ID(p)))

    def processMCParticle(self, p, tuple, prefix, full=False):
        """
        Dump MC Particle to a tuple
        """
        tuple.column(prefix + "TRUE_E", MCE(p))
        tuple.column(prefix + "TRUE_TAU", MCCTAU(p))
        tuple.column(prefix + "TRUE_M", MCM(p) )
        midfun = MCMOTHER( MCID , -1e6 ) 
        gdmidfun = MCMOTHER ( MCMOTHER( MCID , -1e6 ), -1e6)
        motherid = midfun(p)
        gdmotherid = gdmidfun(p)

        if motherid != -1e6 and motherid == gdmotherid:
            print "=================> ", motherid, " <=> ", gdmotherid, " <= "
        tuple.column(prefix + "TRUE_ID", MCID(p) )
        tuple.column(prefix + "TRUE_MOTHER_ID", motherid )
        tuple.column(prefix + "TRUE_GD_MOTHER_ID", gdmotherid )

    def processPartList(self, l, tuple, prefix, NBTRACKS, mcd0top):
        """
        Dump particle list to a tuple
        """
        # Maps with particles properties to keep particle properties
        vctdouble     = cpp.vector('double')
        var_names = [ 'pid', 'pt', 'px', 'py', 'pz', 
                      'e', 'eta', 'phi', 'pt', 'M', 'R', 'vx', 'vy', 'vz',
                      'bpvVD', 'bpvVDChi2', 'bpvIPChi2', 'PIDK',
                      'PIDp', 'PIDpi', 'PIDmu', 'MCID', 'MCMOTHERID', 'MCGMOTHERID', 'MCGGMOTHERID',
                      'D0SAMETREE', 'FSX', 'FSY', 'FSZ']
        
        var_dict = {}
        for vn in var_names:
            var_dict[vn] = vctdouble()

        for (index, p) in enumerate(l):
            var_dict['pid'].push_back(int (ID(p)))
            var_dict['px'].push_back(PX(p))
            var_dict['py'].push_back(PY(p))
            var_dict['pz'].push_back(PZ(p))
            var_dict['e'].push_back(E(p))
            var_dict['eta'].push_back(ETA(p))
            var_dict['phi'].push_back(PHI(p))
            var_dict['pt'].push_back(PT(p))
            var_dict['M'].push_back(MM(p))
            var_dict['vx'].push_back(VX(p.endVertex()))
            var_dict['vy'].push_back(VY(p.endVertex()))
            var_dict['vz'].push_back(VZ(p.endVertex()))
            var_dict['bpvVD'].push_back(BPVVD(p))
            var_dict['bpvVDChi2'].push_back(BPVVDCHI2(p))
            bpvIPCHI2 = BPVIPCHI2 ("") 
            var_dict['bpvIPChi2'].push_back(bpvIPCHI2(p))
            var_dict['PIDK'].push_back(PIDK(p))
            var_dict['PIDp'].push_back(PIDp(p))
            var_dict['PIDpi'].push_back(PIDpi(p))
            var_dict['PIDmu'].push_back(PIDmu(p))
            mcpart = self._getRelatedMCParticle(p)
            midfun = MCMOTHER( MCID , -1e6 )
            gdmidfun = MCMOTHER ( MCMOTHER( MCID , -1e6 ), -1e6)
            ggdmidfun = MCMOTHER ( MCMOTHER ( MCMOTHER( MCID , -1e6 ), -1e6), 1e-6)
            var_dict["MCID"].push_back(MCID(mcpart) )
            var_dict["MCMOTHERID"].push_back(midfun(mcpart) )
            var_dict["MCGMOTHERID"].push_back(gdmidfun(mcpart))
            var_dict["MCGGMOTHERID"].push_back(gdmidfun(mcpart))
            mctop = self._findMCTop(mcpart)
            d0sametree = 0
            if mctop != None and mcd0top != None and mctop == mcd0top:
                d0sametree = 1	
            #print "===> %s D0SameTree: %s" % (prefix, d0sametree)
            var_dict["D0SAMETREE"].push_back(d0sametree)
            fsx = -1e6
            fsy = -1e6
            fsz = -1e6
            if p.proto() != None:
                if p.proto().track() != None:
                    fs = p.proto().track().firstState()
                    if fs != None:
                        fsx = fs.x()
                        fsy = fs.y()
                        fsz = fs.z()
            var_dict['FSX'].push_back(fsx)
            var_dict['FSY'].push_back(fsy)
            var_dict['FSZ'].push_back(fsz)

        # Now add the values to the ntuple
        for k in var_names:
            tuple.farray ( prefix + k , var_dict[k] , prefix + 'nParts' , NBTRACKS)

    def processD0(self, d0, tuple):
        """
        Process a D0 candidate.
        """
        ## Looking for MC particle
        mcd0 = self._getRelatedMCParticle(d0)
        mcd0top = self._findMCTop(mcd0)

        # Just a check
        if mcd0top == None and mcd0 != None:
            print "ERROR ERROR ERROR: This shouldn't happen !"
        # Ignoring ghosts...
        if mcd0top == None or mcd0 == None:
            return

        ## Adding the D0 info itself
        self.processParticle(d0, tuple, "D_", True, firstState=True)

        # Processing the MC
        tuple.column("D_FROMB",  self._isFromB(mcd0))
        self.processMCParticle(mcd0, tuple, 'D_')

        # Checking D0 lifetime        
        d0pv = self.getRelatedPV(d0)
        myctau = CTAU (d0pv)
        tuple.column("D_CTAU",  myctau(d0))

        # Processing daugthers
        for i, d in enumerate(d0.daughters()):
            #name = d.name().replace("+", "").replace("-", "")
            name = "DChild%d" % (i + 1)
            self.processParticle(d, tuple, "%s_" % name, firstState=True)
            self.processParticlePID(d, tuple, "%s_" % name)
            mcpart = self._getRelatedMCParticle(d0)
            self.processMCParticle(mcpart, tuple, '%s_' % name)

        # Now listing tracks
        daughterTracks = [ t.proto().track() for t in d0.daughters() if t.proto() != None ]
        partsFromPV = [ p for p  in self.get('Phys/StdAllNoPIDsPions/Particles')
                        if  self.getRelatedPV(p) == d0pv \
                        and (p.proto() != None and p.proto().track() not in daughterTracks) ]

        partsFromPVSameDecay = []
        partsFromPVOtherDecay = []
        for p in partsFromPV:
            mcp = self._getRelatedMCParticle(p)
            mcptop = self._findMCTop(mcp)
            if mcp == None or mcptop == None:
                continue
            if mcptop == mcd0top:
                partsFromPVSameDecay.append(p)
            else:
                partsFromPVOtherDecay.append(p)
                
        # Sort the particles py pt
        # Take the top PT particles and keep them
        NBTRACKS = 10
        ptsorts = sorted(partsFromPVSameDecay, key=lambda p: p.pt(), reverse=True)
        self.processPartList(ptsorts[:NBTRACKS], tuple, "hpt_same_", NBTRACKS, mcd0top)

        ptsorto = sorted(partsFromPVOtherDecay, key=lambda p: p.pt(), reverse=True)
        self.processPartList(ptsorto[:NBTRACKS], tuple, "hpt_other_", NBTRACKS, mcd0top)

        # Now finding the most displaced vertices
        distsorts = sorted(partsFromPVSameDecay, key=lambda p: BPVVDCHI2(p), reverse=True)
        self.processPartList(distsorts[:NBTRACKS], tuple, "vdchi2_same_", NBTRACKS, mcd0top)

        distsorto = sorted(partsFromPVOtherDecay, key=lambda p: BPVVDCHI2(p), reverse=True)
        self.processPartList(distsorto[:NBTRACKS], tuple, "vdchi2_other_", NBTRACKS, mcd0top)


    ## the main 'analysis' method 
    ############################################################################
    def analyse( self ) :   ## IMPORTANT! 
        """
        The main 'analysis' method
        """        
        ## use the native Python print to stdout:
        #self.Print( 'Now filtering the tracks')
        tuple = self.nTuple('DecayTree')

        ## Find the matching D0
        selection =  self.get('/Event/Phys/SelD2KPi/Particles')
        if selection == None:
            return SUCCESS

        ## Looking up event data
        rc_summary   = self.get('/Event/Rec/Summary').summaryData()
        odin  = self.get_ ( '/Event/DAQ/ODIN'   , True )
        
        ## Process D0 in the event
        for recd0 in selection:

            # First checking whether we have a ghost
            mcd0 = self._getRelatedMCParticle(recd0)
            if mcd0 == None:
                #print "Aborting as d0 is ghost"
                continue

            #print "D0 is true MC"
            # Now dumping d0 info
            self.processD0(recd0, tuple)

            ## put rec-summary information into n-tuple
            self.addRecSummary(tuple, rc_summary)
            ## get ODIN from TES
            tuple.column_aux ( odin ) 


            ## Finally write the data to the tuple
            tuple.write()
            
        return SUCCESS      ## IMPORTANT!!!

    
# =============================================================================

# =============================================================================
## The configuration of the job
def configure ( inputdata        ,    ## the list of input files  
                catalogs = []    ,    ## xml-catalogs (filled by GRID)
                castor   = False ,    ## use the direct access to castor/EOS ? 
                params   = {}    ) :

    ## configure  Track <--> MC relation table  
    import LoKiPhysMC.Track2MC_Configuration
    import LoKiMC.MC
    
    ## import DaVinci 
    from Configurables import DaVinci, GaudiSequencer
    ## delegate the actual configurtaion to DaVinci 
    dv = DaVinci ( DataType   = '2011' ,
                   InputType  = 'MDST',
                   Lumi = True,
                   Simulation = True,
                   DDDBtag="MC11-20111102",
                   CondDBtag="sim-20111111-vc-md100",
                   HistogramFile = "mcd02kpi_tracks7_histo.root",
                   TupleFile = "mcd02kpi_tracks7_ntuple.root",
                   PrintFreq = 1000)
    

    from Configurables import DecayTreeTuple, FilterDesktop, TupleToolGeometry, CombineParticles
    from Configurables import MCDecayTreeTuple, TupleToolMCTruth, MCTupleToolHierarchy
    from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence, DataOnDemand
    from Configurables import CheckPV
    

    # First using CombineParticle to create the D0
    ################################################################################
    #from StandardParticles import  StdAllNoPIDsPions, StdAllNoPIDsKaons
    _pions = DataOnDemand(Location='Phys/StdAllNoPIDsPions/Particles')
    _kaons = DataOnDemand(Location='Phys/StdAllNoPIDsKaons/Particles')

    _d2kpi = CombineParticles("d2kpi")
    _d2kpi.DecayDescriptor = "[D0 -> K- pi+]cc"
    _d2kpi.DaughtersCuts = { "K-"  : "(PT > 500.0) & (0.0 < PIDK)",
                             "pi+" : "(PT > 500.0) & (5.0 > PIDK)",
                             "K+"  : "(PT > 500.0) & (0.0 < PIDK)",
                             "pi-" : "(PT > 500.0) & (5.0 > PIDK) " }
    _d2kpi.MotherCut = "(VFASPF(VCHI2/VDOF)<10)"
    _d2kpi.CombinationCut = "(ADAMASS('D0') < 50.0)"
    _d2kpi.Preambulo = [ 
        "from LoKiPhysMC.decorators import *" ,
        "from PartProp.Nodes import CC"      ]
    #_d2kpi.ReFitPVs = True

    SelD2KPi = Selection( "SelD2KPi",
                          Algorithm= _d2kpi,
                          RequiredSelections=[_pions,_kaons] ) 
    
    SeqD2KPi = SelectionSequence('SeqD2KPi',TopSelection = SelD2KPi)

    
    # Now the CheckPV method to filter algorithms
    c = CheckPV("OnePV")
    c.MinPVs = 1

    # And a sequencer to put them together
    gseq = GaudiSequencer()
    gseq.Members = [ c, SeqD2KPi.sequence() ]
    
    ## define the input data
    setData  ( inputdata , catalogs , castor )
    
    ## get/create application manager
    gaudi = appMgr() 
    
    #
    ## modify/update the configuration:
    #
    
    ## (1) create the algorithm
    alg = TrackFilter( 'TrackFilter' )
    
    #seq = createSequencer()
    ## (2) replace the list of top level algorithm by
    #     new list, which contains only *THIS* algorithm
    gaudi.setAlgorithms( [ gseq, alg ] )
             
    return SUCCESS




# =============================================================================

# =============================================================================
## Job steering 
if __name__ == '__main__' :

    logger.info ( 80*'*'  ) 
    logger.info ( __doc__ ) 
    logger.info ( ' Author  : %s ' %  __author__  ) 
    logger.info ( ' Version : %s ' %  __version__ ) 
    logger.info ( ' Date    : %s ' %  __date__    ) 
    logger.info ( 80*'*'  ) 

    ## job configuration
    ## BKQuery('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/27163003/ALLSTREAMS.DST')
    # inputdata = [
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000001_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000002_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000003_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000004_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000005_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000006_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000007_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000008_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000009_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000010_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000011_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000012_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000013_1.allstreams.dst',
    #     '/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000014_1.allstreams.dst',
    #     ]
  
    inputdata = [ 
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000059_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000028_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000025_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000075_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000031_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000001_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000198_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000121_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000138_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000048_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000087_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000065_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000079_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000190_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000048_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000123_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00014710/0000/00014710_00000193_1.allstreams.dst',
      'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/MC11a/ALLSTREAMS.DST/00015514/0000/00015514_00000075_1.allstreams.dst',
      ]



 
    configure( inputdata , castor = True )
    
    ## event loop
    run(5000)
    #run(5000)
        
# =============================================================================
# The END
# =============================================================================


