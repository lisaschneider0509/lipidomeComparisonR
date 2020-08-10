### R code from vignette source 'ALLintro.Rnw'

###################################################
### code chunk number 1: ALLintro.Rnw:56-59
###################################################
library(ALL)
data(ALL)
show(ALL)


###################################################
### code chunk number 2: ALLintro.Rnw:62-63
###################################################
print(summary(pData(ALL)))


###################################################
### code chunk number 3: ALLintro.Rnw:64-65
###################################################
hist(cvv <- apply(exprs(ALL),1,function(x)sd(x)/mean(x)))


###################################################
### code chunk number 4: ALLintro.Rnw:66-71
###################################################
ok <- cvv > .08 & cvv < .18
fALL <- ALL[ok,]
show(fALL)

allx2 <- data.frame(t(exprs(fALL)), class=ALL$BT)


###################################################
### code chunk number 5: ALLintro.Rnw:73-77
###################################################
library(rpart)
rp1 <- rpart(class~.,data=allx2)
plot(rp1)
text(rp1)


