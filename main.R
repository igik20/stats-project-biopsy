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
biopsy <- biopsy[,-1]

names(biopsy) <- c("ID", "Thickness", "SizeUniformity", "ShapeUniformity", "Adhesion", "CellSize", "Nuclei", "Chromatin", "Nucleoli", "Mitoses", "Verdict")
biopsy$Verdict[biopsy$Verdict == "benign"] <- 0
biopsy$Verdict[biopsy$Verdict == "malignant"] <- 1
biopsy$Verdict <- as.numeric(biopsy$Verdict)
biopsy=drop_na(biopsy)




# Exploratory plots
plot(biopsy[,2:10])

# Predictor correlation analysis
attach(biopsy)
test1 <- cor.test(SizeUniformity, ShapeUniformity)
test1$p.value

corrs <- rcorr(as.matrix(biopsy[,2:10]))
corrs
corrplot(corrs$r, order="hclust", tl.col="firebrick2", tl.srt = 45)

plot(jitter(biopsy$Thickness), jitter(biopsy$Verdict), col=c('#8eca74', '#ffc23f')[as.factor(biopsy$Verdict)], main="Thickness vs. Verdict", xlab="Thickness", ylab="Verdict", pch=1)
legend("topleft", legend=c("benign", "malignant"),
       col=c("#8eca74", "#ffc23f"), pch=1, cex=0.8)




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
for(pred in names(biopsy)[2:6]){
  model <- make.model.1(biopsy$Verdict, pred)
  delta <- model$null.deviance - model$deviance
  deltas %>% add_row(predictor = pred, delta = delta)
}

# Single predictor models
model.thickness <- glm(Verdict~Thickness, family="binomial", data=biopsy)
summary(model.thickness) # ResDev 458.5
model.sizeunif <- glm(Verdict~SizeUniformity, family="binomial", data=biopsy)
summary(model.sizeunif) # ResDev 254.8
model.shape <- glm(Verdict~ShapeUniformity, family="binomial", data=biopsy)
summary(model.shape) # ResDev 267.6
model.adhesion <- glm(Verdict~Adhesion, family="binomial", data=biopsy)
summary(model.adhesion) # ResDev 463.3
model.cellsize <- glm(Verdict~CellSize, family="binomial", data=biopsy)
summary(model.cellsize) # ResDev 452.9
model.nuclei <- glm(Verdict~Nuclei, family="binomial", data=biopsy)
summary(model.nuclei) # ResDev 340.6
model.chromatin <- glm(Verdict~Chromatin, family="binomial", data=biopsy)
summary(model.chromatin) # ResDev 388.2
model.nucleoli <- glm(Verdict~Nucleoli, family="binomial", data=biopsy)
summary(model.nucleoli) # ResDev 464.3
model.mitoses <- glm(Verdict~Mitoses, family="binomial", data=biopsy)
summary(model.mitoses) # ResDev 717.52

single.models <- data.frame(pred=names(biopsy)[2:10],resdev=c(458.5, 254.8, 267.6, 463.3, 452.9, 340.6, 388.2, 464.3, 717.52))
barplot(single.models$resdev, names.arg = single.models$pred, xlab="Predictor", ylab="Residual devaince", main="Single variable predictors", col="lightcoral")

# Manual build up
model1 = glm(Verdict~Thickness,family="binomial",data=biopsy)
summary(model1)
model2 = glm(Verdict~Thickness+SizeUniformity,family="binomial",data=biopsy)
summary(model2)
model3 = glm(Verdict~Thickness+SizeUniformity+ShapeUniformity,family="binomial",data=biopsy)
summary(model3)
model4 = glm(Verdict~Thickness+SizeUniformity*ShapeUniformity,family="binomial",data=biopsy)
summary(model4)
model5 = glm(Verdict~.,family="binomial",data=biopsy)
summary(model5)

coef(model5)

list(model1, model2, model3, model4, model5)
anova(model1, model2, model3, model4, model5, test = "Chisq")
# By now, model4 is the best based on deviance as it has the lowest compared to the other models.
