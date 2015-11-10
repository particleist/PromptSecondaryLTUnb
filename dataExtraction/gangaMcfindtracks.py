#
# Exctraction script for MC d02Kpi
#

script   = 'mcfindtracksFinal.py'

t = JobTemplate( application = Bender( version = "v25r3", module = script ))
t.name = "MCFT07"

j = Job( t, backend = Dirac() )
#j = Job( t, backend = Local() )

bkquery = BKQuery('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-EmNoCuts/Sim05/Trig0x40760037Flagged/Reco12a/Stripping17NoPrescalingFlagged/27163003/ALLSTREAMS.DST')
datafiles = bkquery.getDataset()
j.inputdata = datafiles

j.splitter=SplitByFiles(filesPerJob=10, ignoremissing=False, maxFiles=None)
NTUPLE_NAME = "mcd02kpi_tracks7_ntuple.root"
ntuple = DiracFile(NTUPLE_NAME)
histo = DiracFile("mcd02kpi_tracks7_histo.root")
j.outputfiles = [ntuple, histo] 

j.submit()
