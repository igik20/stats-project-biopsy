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
biopsy <- drop_na(biopsy)

# Exploratory plots
plot(biopsy[,2:10])
boxplot(biopsy[,2:10])
for(var in names(biopsy)[2:10])
{
  hist(biopsy[,var], main=var, col="skyblue")
}

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
summary(model.mitoses) # ResDev 717.5

single.models <- data.frame(pred=names(biopsy)[2:10],resdev=c(458.5, 254.8, 267.6, 463.3, 452.9, 340.6, 388.2, 464.3, 717.5))
barplot(single.models$resdev, names.arg = single.models$pred, xlab="Predictor", ylab="Residual devaince", main="Single variable predictors", col="lightcoral")

# Two-variable models
model2.thickness <- glm(Verdict~SizeUniformity+Thickness, family="binomial", data=biopsy)
summary(model2.thickness) # ResDev 196.6
model2.shape <- glm(Verdict~SizeUniformity+ShapeUniformity, family="binomial", data=biopsy)
summary(model2.shape) # ResDev 221.1
model2.adhesion <- glm(Verdict~SizeUniformity+Adhesion, family="binomial", data=biopsy)
summary(model2.adhesion) # ResDev 229.0
model2.cellsize <- glm(Verdict~SizeUniformity+CellSize, family="binomial", data=biopsy)
summary(model2.cellsize) # ResDev 239.5
model2.nuclei <- glm(Verdict~SizeUniformity+Nuclei, family="binomial", data=biopsy)
summary(model2.nuclei) # ResDev 166.3
model2.chromatin <- glm(Verdict~SizeUniformity+Chromatin, family="binomial", data=biopsy)
summary(model2.chromatin) # ResDev 207.3
model2.nucleoli <- glm(Verdict~SizeUniformity+Nucleoli, family="binomial", data=biopsy)
summary(model2.nucleoli) # ResDev 223.0
model2.mitoses <- glm(Verdict~SizeUniformity+Mitoses, family="binomial", data=biopsy)
summary(model2.mitoses) # ResDev 241.4

double.models <- data.frame(pred2=c("Thickness", "ShapeUniformity", "Adhesion", "CellSize", "Nuclei", "Chromatin", "Nucleoli", "Mitoses"), resdev=c(196.6, 221.1, 229.0, 239.5, 166.3, 207.3, 223.0, 241.4))
barplot(double.models$resdev, names.arg=double.models$pred2, xlab="Second predictor", ylab="Residual deviance", main="Second predictors", col="olivedrab2")

model3.thickness <- glm(Verdict~SizeUniformity+Nuclei+Thickness, family="binomial", data=biopsy)
summary(model3.thickness) # ResDev 135.6
model3.shape <- glm(Verdict~SizeUniformity+Nuclei+ShapeUniformity, family="binomial", data=biopsy)
summary(model3.shape) # ResDev 153.1
model3.adhesion <- glm(Verdict~SizeUniformity+Nuclei+Adhesion, family="binomial", data=biopsy)
summary(model3.adhesion) # ResDev 161.0
model3.cellsize <- glm(Verdict~SizeUniformity+Nuclei+CellSize, family="binomial", data=biopsy)
summary(model3.cellsize) # ResDev 161.6
model3.chromatin <- glm(Verdict~SizeUniformity+Nuclei+Chromatin, family="binomial", data=biopsy)
summary(model3.chromatin) # ResDev 151.8
model3.nucleoli <- glm(Verdict~SizeUniformity+Nuclei+Nucleoli, family="binomial", data=biopsy)
summary(model3.nucleoli) # ResDev 150.9
model3.mitoses <- glm(Verdict~SizeUniformity+Nuclei+Mitoses, family="binomial", data=biopsy)
summary(model3.mitoses) # ResDev 158.5

triple.models <- data.frame(pred3=c("Thickness", "ShapeUniformity", "Adhesion", "CellSize", "Chromatin", "Nucleoli", "Mitoses"), resdev=c(135.6, 153.1, 161.0, 161.6, 151.8, 150.9, 158.5))
barplot(triple.models$resdev, names.arg=triple.models$pred3, xlab="Third predictor", ylab="Residual deviance", main="Third predictors", col="royalblue1")

model4.shape <- glm(Verdict~SizeUniformity+Nuclei+Thickness+ShapeUniformity, family="binomial", data=biopsy)
summary(model4.shape) # ResDev 129.5
model4.adhesion <- glm(Verdict~SizeUniformity+Nuclei+Thickness+Adhesion, family="binomial", data=biopsy)
summary(model4.adhesion) # ResDev 126.3
model4.cellsize <- glm(Verdict~SizeUniformity+Nuclei+Thickness+CellSize, family="binomial", data=biopsy)
summary(model4.cellsize) # ResDev 131.5
model4.chromatin <- glm(Verdict~SizeUniformity+Nuclei+Thickness+Chromatin, family="binomial", data=biopsy)
summary(model4.chromatin) # ResDev 123.3
model4.nucleoli <- glm(Verdict~SizeUniformity+Nuclei+Thickness+Nucleoli, family="binomial", data=biopsy)
summary(model4.nucleoli) # ResDev 123.8
model4.mitoses <- glm(Verdict~SizeUniformity+Nuclei+Thickness+Mitoses, family="binomial", data=biopsy)
summary(model4.mitoses) # ResDev 133.0

quad.models <- data.frame(pred4 = c("ShapeUniformity", "Adhesion", "CellSize", "Chromatin", "Nucleoli", "Mitoses"), resdev=c(129.5, 126.3, 131.5, 123.3, 123.8, 133.0))
barplot(quad.models$resdev, names.arg=quad.models$pred4, xlab="Fourth predictor", ylab="Residual devaince", main="Fourth predictors", col="goldenrod1")

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
