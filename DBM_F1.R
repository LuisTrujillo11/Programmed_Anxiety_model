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

dataPCA <- read.csv("Doctorado2020/Resultados_2020/PCA/PCAF1F2/PCA_database_F1.csv")
dataPCA$jacobians <- paste ("Doctorado2020/Resultados_2020/MRI/F1_preproc/Proc_1levl/jacobians/overall/secondlevel_sub-0",dataPCA$id,"_FLASH_pp_relative.mnc", sep = "")
dataPCA <- dataPCA  %>% slice(1, 3, 4, 10, 11, 13, 14, 18, 20, 22, 24, 25, 30, 31 )
dataPCA$group  <- factor(dataPCA$group)
dataPCA$group <- relevel(dataPCA$group, ref = "CON-NA")
dataPCA$weight <- dataF1$W2mths
# Mask and atlas ----------------------------------------------------------

mask <- "Data/fisher atlas/Fischer344_nifti/Fischer344_mask.mnc"
anat_volume <- mincGetVolume("Data/fisher atlas/Fischer344_nifti/Fischer344_template.mnc")
anat_volume2 <- mincArray(mincGetVolume("Doctorado2020/Resultados_2020/MRI/F1_preproc/Proc_1levl/secondlevel_template0.mnc"))
labels= "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv"
atlas = mincArray(anat_volume)
alaslabel = "Data/fisher atlas/Fischer344_nifti/Fischer344_labels.mnc"


# General lineal model contrast vs CON-NA ----------------------------------------------------

vs <- mincLm(jacobians ~ group + W2mths, dataPCA, mask = mask)
vs_corrdiet_r <- mincFDR(vs, mask = mask)
peaksCON_A <- mincFindPeaks(vs, "tvalue-groupCON-A", minDistance = 5, threshold = 8.827107)
peaksCAF_NA <- mincFindPeaks(vs, "tvalue-groupCAF-NA", minDistance = 5, threshold = 7.352325)
peaksCAF_A <- mincFindPeaks(vs, "tvalue-groupCAF-A", minDistance = 5, threshold = 6.549562)
peaksCON_A <- mincLabelPeaks(peaksCON_A, alaslabel, defs = "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv")
peaksCAF_NA <- mincLabelPeaks(peaksCAF_NA, alaslabel, defs = "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv")
peaksCAF_A <- mincLabelPeaks(peaksCAF_A, alaslabel, defs = "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv")
launch_shinyRMINC(vs, anatVol = anat_volume2, fdr = vs_corrdiet_r, anatLow = 1, anatHigh = 10)

#save significative peaks and tmaps for each contrast

write.csv(peaksCON_A, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/peaksCON_NA_vs_CON_A.csv")
write.csv(peaksCAF_NA, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/peaksCON_NA_vs_CAF_NA.csv")
write.csv(peaksCAF_A, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/peaksCON_NA_vs_CAF_A.csv")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/tmap_F1_CON_NA_vs_CON_A.mnc", column = "tvalue-groupCON-A")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/tmap_F1_CON_NA_vs_CAF_NA.mnc", column = "tvalue-groupCAF-NA")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/tmap_F1_CON_NA_vs_CAF_A.mnc", column = "tvalue-groupCAF-A")

# General lineal model contrast vs CON-A ----------------------------------------------------

dataPCA$group <- relevel(dataPCA$group, ref = "CON-A")
vs <- mincLm(jacobians ~ group, dataPCA, mask = mask)
vs_corrdiet_r <- mincFDR(vs, mask = mask)
launch_shinyRMINC(vs, anatVol = anat_volume2, fdr = vs_corrdiet_r, anatLow = 1, anatHigh = 10, plotcolumns = dataPCA)
peaksCAF_A <- mincFindPeaks(vs, "tvalue-groupCAF-A", minDistance = 5, threshold = 6.482221)
peaksCAF_A <- mincLabelPeaks(peaksCAF_A, alaslabel, defs = "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv")
write.csv(peaksCAF_A, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/peaksCON_A_vs_CAF_A.csv")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F1_preproc/tmap_F1_CON_A_vs_CAF_A.mnc", column = "tvalue-groupCAF-A")

peaksCAF_NA <- mincFindPeaks(vs, "tvalue-groupCAF-NA", minDistance = 5, threshold = 12.380282)
peaksCAF_NA <- mincLabelPeaks(peaksCAF_NA, alaslabel, defs = "Data/fisher atlas/Fischer344_nifti/Fischer344_Label_Hierarchy2.csv")
write.csv(peaksCAF_NA, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/peaksCON_A_vs_CAF_NA.csv")
mincWriteVolume(vs, "Doctorado2020/Resultados_2020/MRI/F1_preproc/mri_phenotype/tmap_F1_CON_A_vs_CAF_NA.mnc", column = "tvalue-groupCAF-NA")


# Generate high resolution MRI images of specific peaks -------------------



poscolours = colorRampPalette(c("red", "yellow"))(255)
negcolours = colorRampPalette(c("blue", "turquoise1"))(255)
graycolours = gray.colors(255)

opar <- par(bg=graycolours[1]) # set the background to be the same colour as the under colour

##CON NA vs CON A F1 

mincImage(mincArray(anatVol), slice=220, axes=F, low=1.4, high=3) #right piriform cortex
mincImage(mincArray(vs, "tvalue-groupCON-A"), slice=220, low=6, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCON-A"), slice=220, low=-5, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-5"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("5", "13"), poscolours, col="white", align="rb")

##CON NA vs CAF A F1 
mincImage(mincArray(anatVol), slice=105, axes=F, low=1.4, high=3) #cerebello 3
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=105, low=6.54, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=105, low=-6.54, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-6.54"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("6.54", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=185, axes=F, low=1.4, high=3) #right thalamus
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=185, low=4, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=185, low=-4, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-6.54"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("6.54", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=63, axes=F, low=1.4, high=3) #crus 1 anisoform lobule
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=63, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=185, axes=F, low=1.4, high=3) #right hippocampus CA3
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=185, low=4, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=185, low=-4, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-6.54"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("6.54", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=318, axes=F, low=1.4, high=3) #left frontal association cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=318, levels=c(5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=63, axes=F, low=1.4, high=3) #cerebellar lobule 6
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=63, levels=c(5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=250, axes=F, low=1.4, high=3) #secondary sensory motor FUERA DE CEREBRO
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=250, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

##CON NA vs CAF NA F1 
mincImage(mincArray(anatVol), slice=105, axes=F, low=1.4, high=3) #cerebellar lobule 3
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=105, low=7.35, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=105, low=-7.35, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-7.35"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("7.35", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=375, axes=F, low=1.4, high=3) #right olfactory bulb
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=375, low=6, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=105, low=-6, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-7.35"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("7.35", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=220, axes=F, low=1.4, high=3) #right hypothalamus
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=220, low=3.5, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=220, low=-6, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-7.35"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("7.35", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=335, axes=F, low=1.4, high=3) #left olfactory bulb
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=335, low=3.5, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=335, low=-6, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-7.35"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("7.35", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=318, axes=F, low=1.4, high=3) #left infralimbic cortex
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=318, low=2, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-NA"), slice=318, low=-6, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-7.35"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("7.35", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=190, axes=F, low=1.4, high=3) #left substantia nigra
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=190, levels=c(6), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=73, axes=F, low=1.4, high=3) #right crus 1 ansiform lobule
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=73, levels=c(3.5), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=244, axes=F, low=1.4, high=3) #right thalamus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=244, levels=c(4), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=63, axes=F, low=1.4, high=3) #cerebellar lobule 6
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=63, levels=c(3.5), add=T, lty=2, col="red",  drawlabels=F)


##CON A vs CAF A F1
mincImage(mincArray(anatVol), slice=105, axes=F, low=1.4, high=3) #cerebello 3
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=105, low=6, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=105, low=-6, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-6.48"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("6.48", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=317, axes=F, low=1.4, high=3) #right primary somatosensory cortex
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=317, low=6, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=317, low=-6, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-6.48"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("6.48", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=76, axes=F, low=1.4, high=3) #right crus 1 annisoform lobule
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=76, levels=c(6.482221), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=85, axes=F, low=1.4, high=3) #left paraflocculus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=85, levels=c(6.482221), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=230, axes=F, low=1.4, high=3) #left paraflocculus
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=230, levels=c(6.482221), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=375, axes=F, low=1.4, high=3) #right olfactory bulb
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=375, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=220, axes=F, low=1.4, high=3) #right primary somatosensory cortex 2
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=220, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=220, axes=F, low=1.4, high=3) #right piriform cortex
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=220, low=6, high=13, add=T, col=poscolours, underTransparent = T)
mincImage(mincArray(vs, "tvalue-groupCAF-A"), slice=220, low=-5, high=-13, add=T, col=negcolours, underTransparent = T)
color.legend(20, -20, 130, 0, c("-13", "-6.48"), rev(negcolours), col="white", align="rb")
color.legend(170, -20, 280, 0, c("6.48", "13"), poscolours, col="white", align="rb")

mincImage(mincArray(anatVol), slice=190, axes=F, low=1.4, high=3) #right hippocampus CA3
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=190, levels=c(6.48), add=T, lty=2, col="red", drawlabels=F)

mincImage(mincArray(anatVol), slice=322, axes=F, low=1.4, high=3) #left frontal association cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=322, levels=c(5), add=T, lty=2, col="red", drawlabels=F)

mincImage(mincArray(anatVol), slice=220, axes=F, low=1.4, high=3) #left amygdaloid area
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=220, levels=c(5), add=T, lty=2, col="red", drawlabels=F)

mincImage(mincArray(anatVol), slice=140, axes=F, low=1.4, high=3) #periaqueductal grey
mincContour(abs(mincArray(vs, "tvalue-groupCAF-A")), slice=140, levels=c(5), add=T, lty=2, col="red", drawlabels=F)

##CON A vs CAF NA F1 12.380282

mincImage(mincArray(anatVol), slice=105, axes=F, low=1.4, high=3) #cerebello 3
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=105, levels=c(12.380282), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=375, axes=F, low=1.4, high=3) #right olfactory bulb
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=375, levels=c(3), add=T, lty=2, col="red",  drawlabels=F)

mincImage(mincArray(anatVol), slice=318, axes=F, low=1.4, high=3) #left infralimbic cortex
mincContour(abs(mincArray(vs, "tvalue-groupCAF-NA")), slice=318, levels=c(7), add=T, lty=2, col="red",  drawlabels=F)




