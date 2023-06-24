# Packages -------------------------------------------

setwd("/media/sf_F_DRIVE/")
setwd("/Volumes/RENGOKU/")

library(RMINC)
library(parallel)
library(readxl)
library(dplyr)
library(tidyverse)
library(plotrix)

# Databases ---------------------------------------------------------------

dataPCA <- read.csv("Doctorado2020/Resultados_2020/PCA/PCAF1F2/PCA_database_F2.csv")
dataPCA$jacobians <- paste ("Doctorado2020/Resultados_2020/MRI/F2/BIDS/proc/secondlevel_sub-0",dataPCA$id2,"_FLASH_pp.mnc_relative.mnc", sep = "")
dataPCA <- dataPCA  %>% slice(1:19, 21:22, 25, 26, 28:31, 33, 35:38, 40:49)
dataPCA$group  <- factor(dataPCA$group)
dataPCA$group <- relevel(dataPCA$group, ref = "CON-NA")

# Mask and atlas ----------------------------------------------------------

mask <- "Doctorado2020/Resultados_2020/MRI/F2/secondlevel_otsumask.mnc"
anat_volume <- mincGetVolume("Data/fisher atlas/Fischer344_nifti/Fischer344_template.mnc")
anat_volume2 <- mincArray(mincGetVolume("Doctorado2020/Resultados_2020/MRI/F1_preproc/Proc_1levl/secondlevel_template0.mnc"))
labels= "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv"
atlas = mincArray(anat_volume)
alaslabel = "Data/fisher atlas/Fischer344_nifti/Fischer344_labels.mnc"

# General lineal model contrast vs CON-NA ----------------------------------------------------

vs <- mincLm(jacobians ~ group, dataPCA, mask = mask)
vs_corrdiet_r <- mincFDR(vs, mask = mask)
peaksCON_A <- mincFindPeaks(vs, "tvalue-groupCON-A", minDistance = 5, threshold = 4.121667)
peaksCAF_NA <- mincFindPeaks(vs, "tvalue-groupCAF-NA", minDistance = 5, threshold = 3.845282)
peaksCAF_A <- mincFindPeaks(vs, "tvalue-groupCAF-A", minDistance = 5, threshold = 4.845575)


32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 48 49 50 51 52 53 54 56 58 59 60 61 62 65 66 67 68 69
70 71 72 73 74 75 76 77 78 79 80 81 82

dataPCA <- dataPCA  %>% slice(1:15, 17:23, 25, 27:31, 33:50)











