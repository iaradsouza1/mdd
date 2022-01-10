
library(purrr)
library(IsoformSwitchAnalyzeR)
library(patchwork)

# PLOT DTU TRANSCRIPTS ----------------------------------------------------

# Get all rds files from ISA step 2
files <- list.files("results/ISA", pattern = "pass2", recursive = T, full.names = T)
regions <- sapply(strsplit(files, split = "\\/"), "[[", 3)

# Load files
ls_isa <- map(files, ~ {
  load(.x)
  SwitchList_2
})
names(ls_isa) <- regions


# DTU recurrently altered in each sex -----------------------------------

# SLX1A
slx1a_cg25 <- switchPlotTranscript(
  ls_isa$Cg25, 
  gene = "SLX1A",
  condition1 = "Cg25_CTRL_female",
  condition2 = "Cg25_MDD_female"
)


slx1a_ofc <- switchPlotTranscript(
  ls_isa$OFC, 
  gene = "SLX1A",
  condition1 = "OFC_CTRL_female",
  condition2 = "OFC_MDD_female"
)

# TOX2
tox2_cg25 <- switchPlotTranscript(
  ls_isa$Cg25, 
  gene = "TOX2",
  condition1 = "Cg25_CTRL_female",
  condition2 = "Cg25_MDD_female"
) 


tox2_sub <- switchPlotTranscript(
  ls_isa$Sub, 
  gene = "TOX2",
  condition1 = "Sub_CTRL_female",
  condition2 = "Sub_MDD_female"
) 


p1 <- wrap_plots(slx1a_cg25, slx1a_ofc, ncol = 1, nrow = 2)
ggsave(p1, filename = "results/plots_paper/dtu_female_slx1a.svg", width = 7, height = 5)

p2 <- wrap_plots(tox2_cg25, tox2_sub, ncol = 1, nrow = 2)
ggsave(p2, filename = "results/plots_paper/dtu_female_tox2.svg", width = 7, height = 5)


# DTU in the intersection with GWAS ---------------------------------------

# BTN3A2
btn3a2 <- switchPlotTranscript(
  ls_isa$Cg25, 
  gene = "BTN3A2",
  condition1 = "Cg25_CTRL_female",
  condition2 = "Cg25_MDD_female"
)

# SEZ6
sez6 <- switchPlotTranscript(
  ls_isa$Cg25, 
  gene = "SEZ6",
  condition1 = "Cg25_CTRL_male",
  condition2 = "Cg25_MDD_male"
)

# GNRH
gnrh1 <- switchPlotTranscript(
  ls_isa$dlPFC, 
  gene = "GNRH1",
  condition1 = "dlPFC_CTRL_male",
  condition2 = "dlPFC_MDD_male"
)

# PANX2
panx2 <- switchPlotTranscript(
  ls_isa$Nac, 
  gene = "PANX2",
  condition1 = "Nac_CTRL_male",
  condition2 = "Nac_MDD_male"
)

# CARM1P1
carm1p1 <- switchPlotTranscript(
  ls_isa$OFC, 
  gene = "CARM1P1",
  condition1 = "OFC_CTRL_female",
  condition2 = "OFC_MDD_female"
)

# XRCC3
xrcc3 <- switchPlotTranscript(
  ls_isa$OFC, 
  gene = "XRCC3",
  condition1 = "OFC_CTRL_female",
  condition2 = "OFC_MDD_female"
)




