#Load Pass1 Data
load('results/ISA/objects/pass1/Sub_pass1_0.1.rds')

#Now you should run the external sequence analysis tools with the output fasta files from isoformSwitchAnalysisPart1
#Running second step
SwitchList_2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList_1, 
  dIFcutoff                 = 0.1,  
  n                         = NA,
  removeNoncodinORFs        = FALSE,
  pathToCPATresultFile      = "results/ISA/objects/Sub/sub_cpat.txt",
  pathToNetSurfP2resultFile = "results/ISA/objects/Sub/sub_netsurf.csv",
  pathToPFAMresultFile      = "results/ISA/objects/Sub/sub_pfam.txt",
  pathToSignalPresultFile   = 'results/ISA/objects/Sub/sub_spres_summary.signalp5',
  codingCutoff              = 0.725,
  outputPlots               = T,
  pathToOutput              = "."
)

save(SwitchList_2, file="results/ISA/objects/pass2/Sub_pass2.rds")
