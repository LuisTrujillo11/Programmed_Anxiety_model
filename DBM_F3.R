# Packages -------------------------------------------

setwd("/Volumes/RENGOKU/")

library(RMINC)
library(parallel)
library(readxl)
library(dplyr)
library(tidyverse)
library(plotrix)


# Databases ---------------------------------------------------------------
dataPCA <- read.csv("Doctorado2020/Resultados_2020/f3/database_PCA.csv")
dataPCA$jacobians <- paste ("Doctorado2020/Resultados_2020/MRI/F3/jacobians/secondlevel_sub-",dataPCA$id,"_FLASH_pp_relative.mnc", sep = "")
dataPCA$group  <- factor(dataPCA$group)
dataPCA$group <- relevel(dataPCA$group, ref = "CON-NA")
dataPCA <- dataPCA[-c(15, 20, 24, 29, 30, 36, 40, 43, 45, 51, 53, 55, 56, 57, 66, 69), ] 

# Mask and atlas ----------------------------------------------------------

mask <- "Doctorado2020/Resultados_2020/MRI/F3/template/Fischer_to_maskF3Warped.mnc"
anat_volume2 <- mincArray(mincGetVolume("Doctorado2020/Resultados_2020/MRI/F3/jacobians/secondlevel_template0.mnc_relative.mnc"))
labels= "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv"
alaslabel = "Doctorado2020/Resultados_2020/MRI/F3/template/Fischer_to_templateF3_labelWarped.mnc"

# General lineal model contrast vs CON-NA ----------------------------------------------------

vs <- mincLm(jacobians ~ group, dataPCA, mask = mask)
vs_corrdiet_r <- mincFDR(vs, mask = mask)
peaksCON_A <- mincFindPeaks(vs, "tvalue-groupCON-A", minDistance = 5, threshold = 4.121667)
peaksCAF_NA <- mincFindPeaks(vs, "tvalue-groupCAF-NA", minDistance = 5, threshold = 3.845282)
peaksCAF_A <- mincFindPeaks(vs, "tvalue-groupCAF-A", minDistance = 5, threshold = 4.845575)

#save significative peaks and tmaps for each contrast

write.csv(peaksCON_A, "Doctorado2020/Resultados_2020/MRI/F3/peaksCON_NA_vs_CON_A.csv")
write.csv(peaksCAF_NA, "Doctorado2020/Resultados_2020/MRI/F3/peaksCON_NA_vs_CAF_NA.csv")
write.csv(peaksCAF_A, "Doctorado2020/Resultados_2020/MRI/F3/peaksCON_NA_vs_CAF_A.csv")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F3/tmap_F3_CON_NA_vs_CON_A.mnc", column = "tvalue-groupCON-A")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F3/tmap_F3_CON_NA_vs_CAF_NA.mnc", column = "tvalue-groupCAF-NA")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F3/tmap_F3_CON_NA_vs_CAF_A.mnc", column = "tvalue-groupCAF-A")


# General lineal model contrast vs CON-A ----------------------------------------------------

dataPCA$group <- relevel(dataPCA$group, ref = "CON-A")
vs <- mincLm(jacobians ~ group, dataPCA, mask = mask)
vs_corrdiet_r <- mincFDR(vs, mask = mask)
launch_shinyRMINC(vs, anatVol = anat_volume2, fdr = vs_corrdiet_r, anatLow = 1, anatHigh = 10, plotcolumns = dataPCA)
peaksCAF_A <- mincFindPeaks(vs, "tvalue-groupCAF-A", minDistance = 5, threshold = 4.351709)
write.csv(peaksCAF_A, "Doctorado2020/Resultados_2020/MRI/F3/peaksCON_A_vs_CAF_A.csv")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F3/tmap_F3_CON_A_vs_CAF_A.mnc", column = "tvalue-groupCAF-A")

peaksCAF_NA <- mincFindPeaks(vs, "tvalue-groupCAF-NA", minDistance = 5, threshold = 4.812580)
write.csv(peaksCAF_NA, "Doctorado2020/Resultados_2020/MRI/F3/peaksCON_A_vs_CAF_NA.csv")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F3/tmap_F3_CON_A_vs_CAF_NA.mnc", column = "tvalue-groupCAF-NA")

# Generate high resolution MRI images of specific peaks -------------------

mincImage(mincArray(anat_volume2), slice=136, axes=F, low=1, high=4.5) #third ventricle
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=136, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=338, axes=F, low=1, high=4.5) #copula
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=338, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=265, axes=F, low=1, high=4.5) #rhinal cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=265, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=155, axes=F, low=1, high=4.5) #amygdaloid area
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=155, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=205, axes=F, low=1, high=4.5) #ca3
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=205, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=352, axes=F, low=1, high=4.5) #Cerebellar White Matter / Arbor vita of cerebellum 
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=352, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=185, axes=F, low=1, high=4.5) #Internal capsule
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=185, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=230, axes=F, low=1, high=4.5) #Retrosplenial Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=230, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=185, axes=F, low=1, high=4.5) #lateral ventricle izq
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=185, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=336, axes=F, low=1, high=4.5) #cerebellar lobule 6
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=336, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=89, axes=F, low=1, high=4.5) #Corpus Callosum & associated White Matter
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=89, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)


##CON NA vs CAF NA F3
mincImage(mincArray(anat_volume2), slice=195, axes=F, low=1, high=4.5) #third ventricle
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=195, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=132, axes=F, low=1, high=4.5) #hypothalamus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=132, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=340, axes=F, low=1, high=4.5) #paramedian lobe
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=340, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=230, axes=F, low=1, high=4.5) #right subincular region
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=230, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=309, axes=F, low=1, high=4.5) #flocculus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=309, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=358, axes=F, low=1, high=4.5) #cerebellar lobule 8
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=358, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=65, axes=F, low=1, high=4.5) #Cingulate Cortex area 1
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=65, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=217, axes=F, low=1, high=4.5) #Hippocampal layer CA1
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=217, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=80, axes=F, low=1, high=4.5) #Corpus Callosum & associated White Matter
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=80, levels=c(4.5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=185, axes=F, low=1, high=4.5) #Corpus Callosum & associated White Matter 2
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=185, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=260, axes=F, low=1, high=4.5) #subicular region
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=260, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=185, axes=F, low=1, high=4.5) #Internal capsule
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=185, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=168, axes=F, low=1, high=4.5) #Stria terminalis
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=168, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=253, axes=F, low=1, high=4.5) #Secondary Visual Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=253, levels=c(3.5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=115, axes=F, low=1, high=4.5) #Piriform Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=115, levels=c(2), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=40, axes=F, low=1, high=4.5) #Olfactory Bulb
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=40, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

##CON NA vs CON A F3
mincImage(mincArray(anat_volume2), slice=135, axes=F, low=1, high=4.5) #caudaputamen
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=135, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=180, axes=F, low=1, high=4.5) #caudaputamen 2
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=180, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=71, axes=F, low=1, high=4.5) #Frontal Association Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=71, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=207, axes=F, low=1, high=4.5) #Retrosplenial Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=207, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=157, axes=F, low=1, high=4.5) #amygdaloid area
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=157, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=170, axes=F, low=1, high=4.5) #Thalamus
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=170, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=143, axes=F, low=1, high=4.5) #Primary Somatosensory Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=143, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=320, axes=F, low=1, high=4.5) #Crus 1 ansiform lobule
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=320, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=230, axes=F, low=1, high=4.5) #Corpus Callosum & associated White Matter
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=230, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=350, axes=F, low=1, high=4.5) #Cerebellar White Matter / Arbor vita of cerebellum 
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=350, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=220, axes=F, low=1, high=4.5) #dentate gyrus
mincContour(abs(mincArray(vs, "tvalue-groupCON-A")), slice=220, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

##CON A vs CAF A F3
mincImage(mincArray(anat_volume2), slice=330, axes=F, low=1, high=4.5) #paramedian lobe
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=330, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=140, axes=F, low=1, high=4.5) #Primary Somatosensory Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=140, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=162, axes=F, low=1, high=4.5) #left amygdaloid area
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=162, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=159, axes=F, low=1, high=4.5) #right amygdaloid area 
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=159, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=270, axes=F, low=1, high=4.5) #rhinal cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=270, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=185, axes=F, low=1, high=4.5) #Thalamus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=185, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=75, axes=F, low=1, high=4.5) #Frontal Association Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=75, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=325, axes=F, low=1, high=4.5) #paraflocculus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=325, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=110, axes=F, low=1, high=4.5) #Primary Somatosensory Cortex2
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=110, levels=c(3.5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=160, axes=F, low=1, high=4.5) #left amygdaloid area 2
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=160, levels=c(3.5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=160, axes=F, low=1, high=4.5) #Internal capsule
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=160, levels=c(3.5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=320, axes=F, low=1, high=4.5) #cerebellar lobule 4&5
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=160, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=115, axes=F, low=1, high=4.5) #lateral septum
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=115, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=207, axes=F, low=1, high=4.5) #Retrosplenial Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=207, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=58, axes=F, low=1, high=4.5) #Frontal Association Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=58, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=25, axes=F, low=1, high=4.5) #Frontal Association Cortex2
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=25, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)


##CON A vs CAF NA F3
mincImage(mincArray(anat_volume2), slice=143, axes=F, low=1, high=4.5) #Primary Somatosensory Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=143, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=330, axes=F, low=1, high=4.5) #Cerebellar White Matter
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=330, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=160, axes=F, low=1, high=4.5) #right amygdaloid area
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=160, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=165, axes=F, low=1, high=4.5) #left amygdaloid area
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=165, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=144, axes=F, low=1, high=4.5) #caudaputamen
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=144, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=58, axes=F, low=1, high=4.5) #Frontal Association Cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=58, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anat_volume2), slice=160, axes=F, low=1, high=4.5) #Lateral ventricle
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=160, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)


