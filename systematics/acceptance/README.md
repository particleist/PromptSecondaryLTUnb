D0 to K pi lifetime acceptance checks
=====================================

Check of effect of 4mm radial cut
---------------------------------

The following files are used:

RadialSelector.C, RadialSelector.h: TSelector to iterate over the tree and get the acceptance.
Run.C: Run the TSelector over root://eoslhcb.cern.ch//eos/lhcb/user/m/malexand/d2hh/secondariesTagging/mcd02kpi_tracks9_merged.root
RunProoff.C: Run over the same file using ProofLite
