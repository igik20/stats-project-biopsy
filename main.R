# Biopsy Statistics Project fall 2022
# by Xavier Crespo and Igor Trujnara

# Packages
library("ggplot2")
library("corrplot")
library("Hmisc")
library("tidyr")
library("dplyr")

# Data import
biopsy <- read.csv("biopsy.csv")

# Set human-readable names
names(biopsy) <- c("Number", "ID", "Thickness", "SizeUniformity", "ShapeUniformity", "Adhesion", "CellSize", "Nuclei", "Chromatin", "Nucleoli", "Mitoses", "Verdict")

# Exploratory plots
plot(biopsy[,3:11])

# Predictor correlation analysis
attach(biopsy)
test1 <- cor.test(SizeUniformity, ShapeUniformity)
test1$p.value

corrs <- rcorr(as.matrix(biopsy[,3:11]))
corrs
corrplot(corrs$r, order="hclust", tl.col="firebrick2", tl.srt = 45)
