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

# Set response variable to binary
biopsy$Verdict <- ifelse(biopsy$Verdict=="benign", 0, 1)

# Exploratory plot for thickness
plot(jitter(Thickness), jitter(biopsy$Verdict), col="chartreuse3", main="Thickness vs. Verdict", xlab="Thickness", ylab="Verdict")

# Model for thickness as a sample
model.thickness <- glm(Verdict~Thickness, data=biopsy, family=binomial(link="logit"))
summary(model.thickness)
plot(model.thickness)
curve(predict(model.thickness, data.frame(Thickness = x), type="resp"), add=T, col="chocolate1", lwd=2)

# Function to generate models
make.model.1 <- function(resp, pred) {
  model <- glm(resp~pred, family=binomial(link="logit"))
  nD <- model$null.deviance
  rD <- model$deviance
  dD <- (nD - rD)
  return (model)
}

# Testing the function
model.sunif <- make.model.1(biopsy$Verdict, biopsy$SizeUniformity)
summary(model.sunif)
model.sunif$null.deviance - model.sunif$deviance

deltas <- data.frame("template", 0)
names(deltas) <- c("predictor", "delta")
for(pred in names(biopsy)[3:7]){
  model <- make.model.1(biopsy$Verdict, pred)
  delta <- model$null.deviance - model$deviance
  deltas %>% add_row(predictor = pred, delta = delta)
}
