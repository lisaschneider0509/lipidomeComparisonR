### R code from vignette source 'annaffy.Rnw'

###################################################
### code chunk number 1: annaffy.Rnw:71-72
###################################################
library("annaffy")


###################################################
### code chunk number 2: annaffy.Rnw:78-83
###################################################
data(aafExpr)
probeids <- featureNames(aafExpr)
symbols <- aafSymbol(probeids, "hgu95av2.db")
symbols[[54]]
symbols[55:57]


###################################################
### code chunk number 3: annaffy.Rnw:93-94
###################################################
getText(symbols[54:57])


###################################################
### code chunk number 4: annaffy.Rnw:99-101
###################################################
gos <- aafGO(probeids, "hgu95av2.db")
gos[[3]]


###################################################
### code chunk number 5: annaffy.Rnw:110-111
###################################################
gos[[3]][[1]]@name


###################################################
### code chunk number 6: annaffy.Rnw:142-144
###################################################
gbs <- aafGenBank(probeids, "hgu95av2.db")
getURL(gbs[[1]])


###################################################
### code chunk number 7: annaffy.Rnw:156-158
###################################################
lls <- aafLocusLink(probeids, "hgu95av2.db") 
getURL(lls[[2]])


###################################################
### code chunk number 8: annaffy.Rnw:167-169
###################################################
bands <- aafCytoband(probeids, "hgu95av2.db") 
getURL(bands[[2]])


###################################################
### code chunk number 9: annaffy.Rnw:179-181
###################################################
pmids <- aafPubMed(probeids, "hgu95av2.db")
getURL(pmids[[2]])


###################################################
### code chunk number 10: annaffy.Rnw:190-191
###################################################
getURL(gos[[1]])


###################################################
### code chunk number 11: annaffy.Rnw:198-199
###################################################
getURL(gos[[1]][[4]])


###################################################
### code chunk number 12: annaffy.Rnw:209-211
###################################################
paths <- aafPathway(probeids, "hgu95av2.db")
getURL(paths[[4]])


###################################################
### code chunk number 13: annaffy.Rnw:255-256
###################################################
library(multtest)


###################################################
### code chunk number 14: annaffy.Rnw:264-265
###################################################
class <- as.integer(pData(aafExpr)$covar1) - 1


###################################################
### code chunk number 15: annaffy.Rnw:274-277
###################################################
teststat <- mt.teststat(exprs(aafExpr), class)
index <- order(abs(teststat), decreasing = TRUE)
probeids <- featureNames(aafExpr)[index]


###################################################
### code chunk number 16: annaffy.Rnw:291-292
###################################################
aaf.handler()


###################################################
### code chunk number 17: annaffy.Rnw:298-299
###################################################
anncols <- aaf.handler()[c(1:3,8:9,11:13)]


###################################################
### code chunk number 18: annaffy.Rnw:323-324
###################################################
anntable <- aafTableAnn(probeids[1:50], "hgu95av2.db", anncols)


###################################################
### code chunk number 19: annaffy.Rnw:331-332
###################################################
saveHTML(anntable, "example1.html", title = "Example Table without Data")


###################################################
### code chunk number 20: annaffy.Rnw:346-348
###################################################
testtable <- aafTable("t-statistic" = teststat[index[1:50]], signed = TRUE)
table <- merge(anntable, testtable)


###################################################
### code chunk number 21: annaffy.Rnw:366-369
###################################################
exprtable <- aafTableInt(aafExpr, probeids = probeids[1:50])
table <- merge(table, exprtable)
saveHTML(table, "example2.html", title = "Example Table with Data")


###################################################
### code chunk number 22: annaffy.Rnw:379-380
###################################################
saveText(table, "example2.txt")


###################################################
### code chunk number 23: annaffy.Rnw:406-410
###################################################
library(hgu95av2.db)
probeids <- ls(hgu95av2SYMBOL)
gos <- aafGO(probeids[1:2], "hgu95av2.db")
getText(gos)


###################################################
### code chunk number 24: annaffy.Rnw:418-421
###################################################
kinases <- aafSearchText("hgu95av2.db", "Description", "kinase")
kinases[1:5]
print(length(kinases))


###################################################
### code chunk number 25: annaffy.Rnw:429-432
###################################################
probes <- aafSearchText("hgu95av2.db", c("Description", "Pathway"),
                        c("ribosome", "polymerase"))
print(length(probes))


###################################################
### code chunk number 26: annaffy.Rnw:445-448
###################################################
probes <- aafSearchText("hgu95av2.db", "Description",
                        c("DNA", "polymerase"), logic = "AND")
print(length(probes))


###################################################
### code chunk number 27: annaffy.Rnw:456-458
###################################################
gbs <- c("AF035121", "AL021546", "AJ006123", "AL080082", "AI289489")
aafSearchText("hgu95av2.db", "GenBank", gbs)


###################################################
### code chunk number 28: annaffy.Rnw:482-485
###################################################
aafSearchGO("hgu95av2.db", c("GO:0000002", "GO:0000008"))
aafSearchGO("hgu95av2.db", c("2", "8"))
aafSearchGO("hgu95av2.db", c(2, 8))


