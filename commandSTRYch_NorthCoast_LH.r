###### Chiara Barbieri & Lea Lorene Huber
### June 2025
### commands from the paper Huber et al. 
# Human population history on the North Coast of Peru from Y chromosomes and mitogenomes 


library(plyr)
library(tidyverse)
ydata<-read.table("DatasetContinent17loci.txt", header=T, as.is=T, sep="\t")
# the dataset contains individual haplotypes and information about region and language for each individual
npop <- length(spop) # number of populations
popN <- names(spop) # names of the populations

# mutation rates from YHRD
rates<-read.table("mutationRatesSTRfromYHRD.txt", header=T, as.is=T)

ancientpops = c("Tompullo", "Humahuaca", "ElDiquecito", "Muisca", "Candelaria", "Lauricocha")
coastpops = c("Cao", "Tumbes", "Lambayeque", "Piura")

ca = popTable[which(popTable$perugroup != "OUT"),1] # list with Central Andes populations

# colours
coloursMacroG = c("#1B9E77","#7570B3", "#D95F02" ,"#E7298A", "#66A61E" ,"#E6AB02" ,"#A6761D", "#E41A1C", "#377EB8")
coloursPeruG = c("#1B9E77", "#7570B3","#D95F02","#E7298A", "#E6AB02" ,"#E41A1C")

## dataset with 23 loci
ydata23 = read.table("DatasetContinent23loci.txt", header=T, as.is=T, sep="\t")
ydata23 = filter(ydata23, population != "OUT")
dati23 = select(ydata23, DYS19:YGATAH4)
dati23coast = select(ydata23coast, DYS19:YGATAH4)
ydatainfo = select(ydata, ID, Language_Clustering, ecoregion, perugroup, reference) # add needed columns
ydata23 = plyr::join(ydata23, ydatainfo, by = "ID")
ydata23$population = as.factor(ydata23$population)
#ydata23 = rename(ydata23, Population = population)
ydata23coast = filter(ydata23, population %in% c("Cao", "Lambayeque", "Tumbes", "Piura"))
CoastIDs = filter(ydata23, str_detect(ID, "TUM")|str_detect(ID, "LAM")|str_detect(ID, "PIU")|str_detect(ID, "LAL") )$ID


###################################
# Visualize Samples on a Map

popTable<-read.table("listPopDatasetSTR_2017.txt", sep = "\t", header=T)

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
map.world <- map_data(map="world")
gg <- ggplot()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=c(-47,70), xlim=c(-140,-45))   ## coordinates of lat and lon to frame the map area to plot whole america
gg2<-gg+coord_quickmap(ylim=c(-20,11), xlim=c(-85,-65))  #for central andes region and up to panama

gg2 + geom_label_repel(data=filter(popTable, perugroup != "OUT"),
                       aes(x=lon, y=lat,label=Population, color=MacroGroup), size=1, label.padding=0.1, max.overlaps = 25)+
  geom_point(data=filter(popTable, perugroup != "OUT"), aes(x=lon, y=lat, color=MacroGroup), size=1 ) +
  theme(text = element_text(size=12))+
  scale_color_manual(values = coloursPeruG)+
  theme_pubr(legend = "right")


# computing and saving big distance matrices ------------------------------
#### pairwise distances
## 23 dataset
comp1<-matrix(NA, nrow(ydata23coast), nrow(ydata23coast))
rownames(comp1)<-ydata23coast$ID
colnames(comp1)<-ydata23coast$ID

# fill the matrix with pairwise distances between individual haplotypes. measured as the sum of number of steps of differences between each locus.
# takes a while, about 3 hours
for (i in 1:nrow(ydata23coast)) { 			# columns
  for (j in 1:nrow(ydata23coast)){ 			# rows
    comp1[j,i]<-sum(abs(as.numeric(dati23coast[i,])-as.numeric(dati23coast[j,])),na.rm=TRUE)
  }
}
write.table(comp1,"matrixPairDist23.txt", sep="\t") #makes it with shifted 1st row, just fix it in excel


########################
#### ADJUST FOR MUTATION RATE ON A SET OF 15 LOCI: pairwise TMRCA
# mutation rates from YHRD
rates<-read.table("mutationRatesSTRfromYHRD.txt", header=T, as.is=T)

rates = rates[-16,] # keep only rates for 15 loci

subsetSTR<-ydata1[,which(colnames(ydata1)%in%rates$locus)]
comp3<-matrix(NA, nrow(ydata1), nrow(ydata1))
rownames(comp3)<-ydata1$ID
colnames(comp3)<-ydata1$ID

generationyears<- 30 # the calculations will return the number of generations given the mutation rate, we adjust for number of years assuming that for male humans one generation is equivalent to 30 years

# formula Average Squared Difference
tmrca_asd <- function(hap1, hap2, rates) {
  L <- length(hap1)
  squared_diff <- (hap1 - hap2)^2
  tmrca <- sum(squared_diff / rates) / (2 * L)
}
# Calculate TMRCA (in generations)
tmrca_gen <- tmrca_asd(subsetSTR[i,],subsetSTR[j,], rates$X2022Combined)

# Convert to years (e.g., assuming 30 years per generation)
tmrca_years <- tmrca_gen * generationyears


for (i in 1:nrow(subsetSTR)) { 			#colonne
  for (j in 1:nrow(subsetSTR)){ 			#righe
    comp3[j,i]<-tmrca_asd(subsetSTR[i,],subsetSTR[j,], rates$X2022Combined) * generationyears
 # the formula to calculate the distance between individual sequences in years ago with ASD 
  }
}

### original formula with variance

for (i in 3819:nrow(subsetSTR)) { 			#columns
  for (j in 1:nrow(subsetSTR)){ 			#rows
    varM<-apply(rbind(subsetSTR[i,],subsetSTR[j,]),2,var,na.rm=TRUE) # how does na.rm work, is var good for missing values? 
    comp3[j,i]<-mean(varM/rates$X2022Combined,na.rm=TRUE)*generationyears # the formula to calculate the distance between individual sequences in years ago
  }
}# runs very long, in 2h 1109x4000 individuals
as.matrix(comp3)-> mat11
write.table(mat11,"megamatrixPERU15lociadjustforTMRCApt5.txt" , sep="\t") # makes it with shifted first column but fixed it in excel

# Import of big matrices and computation of pop comparison matrices --------------
### pairwise distances
comp2 = read.table("matrixPairDist1904.txt", sep = "\t", header = TRUE, row.names = 1 )
comp = as.matrix(comp2)


## making matrices of pop distances (needed directly for heatmaps) and saving in melt form
makematavsharing = function(npop, popN, spop, ydata, comp){
  # the average sharing between populations
  matrixaveragesharingPop<-matrix(NA, npop, npop)
  rownames(matrixaveragesharingPop)<-popN
  colnames(matrixaveragesharingPop)<-popN
  # number of identical haplotypes between populations
  matrixNidenticalPop<-matrix(NA, npop, npop)
  rownames(matrixNidenticalPop)<-popN
  colnames(matrixNidenticalPop)<-popN
  # similar haplotypes between populations, arbitrary value less than 2 mutations steps apart
  matrixSimilarPop<-matrix(NA, npop, npop)
  rownames(matrixSimilarPop)<-popN
  colnames(matrixSimilarPop)<-popN
  # similar haplotypes between populations, but with higher distances, arbitrary value less than 4 mutations steps apart
  matrixLESSSimilarPop<-matrix(NA, npop, npop)
  rownames(matrixLESSSimilarPop)<-popN
  colnames(matrixLESSSimilarPop)<-popN
  
  for (i in 1:npop) { 			# columns
    for (j in 1:npop){ 			# rows
      #temp<-comp1[which(rownames(comp1)==popN[j]),which(colnames(comp1)==popN[i])]
      IDpop1<-ydata$ID[which(ydata$population == popN[j])]#take the combination of individuals from pop[i] and pop[j]
      IDpop2<-ydata$ID[which(ydata$population == popN[i])]
      temp<-comp[which(rownames(comp)%in%IDpop1),which(colnames(comp)%in%IDpop2)] #and the distance between the individuals 
      matrixaveragesharingPop[j,i]<-mean(temp)  # mean of the distance values
      matrixNidenticalPop[j,i]<-length(which(temp==0))/(spop[i]*spop[j])   # count how many pairs have distance =0 (identical haplotypes), then adjust for population size
      matrixSimilarPop[j,i]<-length(which(temp<2))/(spop[i]*spop[j])   # adjust for population size
      matrixLESSSimilarPop[j,i]<-length(which(temp<4))/(spop[i]*spop[j])   # adjust for population size
    }
  }
  distlist = list(matrixaveragesharingPop, matrixNidenticalPop, matrixSimilarPop, matrixLESSSimilarPop)
  return(distlist)
}
distlist = makematavsharing(length(spop), names(spop), spop, ydata, comp) # makes different options for ydata possible
matavsharingPop = distlist[[1]]
matrixNidenticalPop = distlist[[2]]
matrixSimilarPop = distlist[[3]]
matrixLESSSimilarPop = distlist[[4]]

## same for 23 loci set: import
matdist23 = read.table("matrixPairDist23.txt", sep = "\t", header = TRUE, row.names = 1)
matdist23 = as.matrix(matdist23)
# sharing between populations
ydata23 = rename(ydata23, population = Population)
distlist23 = makematavsharing(nlevels(ydata23$population),levels(ydata23$population),
                              table(ydata23$population), ydata23, matdist23 )
matavsharing23 = distlist23[[1]]
matNidentical23 = distlist23[[2]]
matSimilar23 = distlist23[[3]]
matLessSimilar23 = distlist23[[4]]

### TMRCA ***********************************
mat<-read.table("megamatrixPERU15lociadjustforTMRCA.txt" , sep="\t", header = TRUE, row.names = 1)

mat<-as.matrix(mat)
mat1 = mat[-which(rownames(mat)%in%arcticpopIDs), -which(colnames(mat)%in%arcticpopIDs)] # exclude Arctic pops
mat1 = mat[-which(rownames(mat)%in%arcticpopIDs), -which(colnames(mat)%in%arcticpopIDs)] # exclude Mocovi and TobaSantaFe
rownames(mat1)<-ydata1$population[match(rownames(mat1),ydata1$ID)]
colnames(mat1)<-ydata1$population[match(colnames(mat1),ydata1$ID)]

# compute between population times etc.
matrixaverageTMRCAPop<-matrix(NA, npop, npop)
rownames(matrixaverageTMRCAPop)<-popN
colnames(matrixaverageTMRCAPop)<-popN

matrixmedianTMRCAPop<-matrix(NA, npop, npop)
rownames(matrixmedianTMRCAPop)<-popN
colnames(matrixmedianTMRCAPop)<-popN

matrixZeroTMRCAPop<-matrix(NA, npop, npop)
rownames(matrixZeroTMRCAPop)<-popN
colnames(matrixZeroTMRCAPop)<-popN

matrix500TMRCA<-matrix(NA, npop, npop)
rownames(matrix500TMRCA)<-popN
colnames(matrix500TMRCA)<-popN

matrix1000TMRCA<-matrix(NA, npop, npop)
rownames(matrix1000TMRCA)<-popN
colnames(matrix1000TMRCA)<-popN

matrix2000TMRCA<-matrix(NA, npop, npop)
rownames(matrix2000TMRCA)<-popN
colnames(matrix2000TMRCA)<-popN

matrix3000TMRCA<-matrix(NA, npop, npop)
rownames(matrix3000TMRCA)<-popN
colnames(matrix3000TMRCA)<-popN

for (i in 1:npop) { 			#COLUMNS
  for (j in 1:npop){ 			#ROWS
    temp<-mat1[which(rownames(mat1)==popN[j]),which(colnames(mat1)==popN[i])] 
    matrixaverageTMRCAPop[j,i]<-mean(temp)
    matrixmedianTMRCAPop[j,i]<-median(temp)
    matrixZeroTMRCAPop[j,i]<-length(which(temp==0))/(spop[i]*spop[j])   # adjust for population size
    matrix500TMRCA[j,i]<-length(which(temp<500))/(spop[i]*spop[j])   # adjust for population size
    matrix1000TMRCA[j,i]<-length(which(temp<1000))/(spop[i]*spop[j])   # adjust for population size
    matrix2000TMRCA[j,i]<-length(which(temp<2000))/(spop[i]*spop[j])   # adjust for population size
    matrix3000TMRCA[j,i]<-length(which(temp<3000))/(spop[i]*spop[j])   # adjust for population size
  }
}

##### ##### ##### ##### ##### ##### 
##### population distances
### import Rst data calculated with YHRD tool ************************************************************
##### ##### ##### ##### ##### ##### 

library(HelpersMG)
matRst = read.table("RstYHRD17loci.txt", sep="\t", header=T)  # calculated with online tool at YHRD
matRst = as.matrix(matRst)
rownames(matRst) = matRst[,1]
matRst= matRst[,-1]
matRst[upper.tri(matRst, diag = TRUE)] = NA
matRst = symmetricize(matRst, method = "ld")
diag(matRst) = 0
class(matRst) = "numeric"
for (i in 1:nrow(matRst)){ # negative Rst values are an artefact of formula, replace with 0 for biological meaning
  for (j in 1:ncol(matRst)){
    if (matRst[i,j] < 0){
      matRst[i,j] = 0
    }
  }
}

# heat maps -----------------------------
library(gplots)
library(RColorBrewer)
library(wesanderson)

heatoranges = c("#FFFFFF", colorRampPalette(brewer.pal(8, "Greens"))(10)) # I want palette to start at pure white
makeheatmap = function(mat, popTable, pdfname, grouping){ # function for making heat maps with dendrograms
  rownames(popTable)<-popTable$Population
  popTable2 = popTable[popTable$Population %in% rownames(mat),]
  popTableSort = data.frame(Population = rownames(mat))
  popSorted = join(popTableSort, popTable2, by = "Population")
  # color the names of the populations according to region criteria
  region<-table(popSorted[[grouping]])
  if(grouping == "MacroGroup"){
    colorcod = coloursMacroG[1:nlevels(popSorted[[grouping]])]
  }else if(grouping == "perugroup"){
    colorcod = coloursPeruG[1:nlevels(popSorted[[grouping]])]
  }
  colorstring<-c(rep(NA, length(popSorted[[grouping]])))
  for (i in 1:nrow(mat)) {
    colorstring[which(popSorted[[grouping]]==labels(region[i]))]<-colorcod[i]
  }
  heatoranges = c("#FFFFFF", colorRampPalette(brewer.pal(8, "Oranges"))(100))# I want palette to start at pure white
  pdf(pdfname)
  h = heatmap.2(mat, scale = "none", col = heatoranges, 
                trace = "none", density.info = "none", cexRow = 0.2, cexCol = 0.2, colRow = colorstring, colCol = colorstring)
  dev.off()
}

## creating different subsets of pop matrices for heat map and generally prepare them********************
# create new grouping for only coast individuals heat map (23 markers)
grouping = "population"
mat = matdist23[(rownames(matdist23)%in%CoastIDs),(rownames(matdist23)%in%CoastIDs)]
rownames(ydata23) = ydata23$ID
ydata23sub = ydata23[ydata23$ID %in% rownames(mat),]
ydata23sub$population = droplevels(ydata23sub$population)
coastNrs = table(ydata23sub[[grouping]]) # population is the grouping here
colorcod = wes_palette(n=nlevels(ydata23sub[[grouping]]), name = "GrandBudapest1")
colorstring = c(rep(NA, length(ydata23sub[[grouping]])))
for (i in 1:nrow(mat)){
  colorstring[which(ydata23sub[[grouping]]==labels(coastNrs[i]))] = colorcod[i]
}
heatoranges = c("#FFFFFF", colorRampPalette(brewer.pal(8, "Oranges"))(100))# I want palette to start at pure white
pdf("heatmapCoastInds.pdf")
h = heatmap.2(mat, scale = "none", col = heatoranges, 
              trace = "none", density.info = "none", cexRow = 0.4, cexCol = 0.4, colRow = colorstring, colCol = colorstring)
dev.off()

# for identical and similar haplotypes, diagonal is distracting so make it white i.e. 0
diag(matrixNidenticalPop) = 0
diag(matrixSimilarPop) = 0
diag(matrix2000TMRCA) = 0
diag(matrix1000TMRCA) = 0
diag(matrix3000TMRCA) = 0
diag(matrix500TMRCA) = 0

matRst1 = matRst[!rownames(matRst)%in%arcticpops,!rownames(matRst)%in%arcticpops]

## generating heat map plots ************************************
makeheatmap(matrix2000TMRCA, popTable, "heatmapTMRCA2000.pdf", "MacroGroup")


# making and saving of meltpop versions -------------------
library(reshape)
### pairwise distances
meltpop3<-melt(matrixaveragesharingPop)  # the same data but in melt form instead of in matrix form
colnames(meltpop3)<-c("pop1", "pop2", "averageDifference")

meltpop3$identical<-melt(matrixNidenticalPop)$value
meltpop3$similar<-melt(matrixSimilarPop)$value
meltpop3$LESSsimilar<-melt(matrixLESSSimilarPop)$value

meltpop3$region1<-ydata$ecoregion[match(meltpop3$pop1,ydata$population)]
meltpop3$region2<-ydata$ecoregion[match(meltpop3$pop2,ydata$population)]
meltpop3$group1<-ydata$perugroup[match(meltpop3$pop1,ydata$population)]
meltpop3$group2<-ydata$perugroup[match(meltpop3$pop2,ydata$population)]

for (i in 1:nrow(meltpop3)){
  ttt<-sort(c(meltpop3$region1[i], meltpop3$region2[i]))
  meltpop3$regionCouple[i]<-paste(ttt, sep="+", collapse="+")
}
for (i in 1:nrow(meltpop3)){
  ttt<-sort(c(meltpop3$group1[i], meltpop3$group2[i]))
  meltpop3$groupCouple[i]<-paste(ttt, sep="+", collapse="+")
}


### TMRCA
meltpop<-melt(matrixaverageTMRCAPop)
colnames(meltpop)<-c("pop1", "pop2", "averageTMRCA")
meltpop$identical<-melt(matrixZeroTMRCAPop)$value
meltpop$TMRCA500<-melt(matrix500TMRCA)$value
meltpop$TMRCA1000<-melt(matrix1000TMRCA)$value
meltpop$TMRCA2000<-melt(matrix2000TMRCA)$value
meltpop$TMRCA3000<-melt(matrix3000TMRCA)$value

meltpop3 = meltpop

meltpop3$region1<-ydata1$ecoregion[match(meltpop3$pop1,ydata1$population)]
meltpop3$region2<-ydata1$ecoregion[match(meltpop3$pop2,ydata1$population)]
meltpop3$group1<-ydata1$perugroup[match(meltpop3$pop1,ydata1$population)]
meltpop3$group2<-ydata1$perugroup[match(meltpop3$pop2,ydata1$population)]
meltpop3$lang1<-ydata1$Language_Clustering[match(meltpop3$pop1,ydata1$population)]
meltpop3$lang2<-ydata1$Language_Clustering[match(meltpop3$pop2,ydata1$population)]

for (i in 1:nrow(meltpop3)){
  ttt<-sort(c(meltpop3$region1[i], meltpop3$region2[i]))
  meltpop3$regionCouple[i]<-paste(ttt, sep="+", collapse="+")
}
for (i in 1:nrow(meltpop3)){
  ttt<-sort(c(meltpop3$group1[i], meltpop3$group2[i]))
  meltpop3$groupCouple[i]<-paste(ttt, sep="+", collapse="+")
}
for (i in 1:nrow(meltpop3)){
  ttt<-sort(c(meltpop3$lang1[i], meltpop3$lang2[i]))
  meltpop3$langCouple[i]<-paste(ttt, sep="+", collapse="+")
}

### Rst meltpop for boxplots (saved them, import below) don't run
meltpopRst<-melt(matRst)
colnames(meltpopRst)<-c("pop1", "pop2", "Rst")

meltpopRst$region1<-ydata$ecoregion[match(meltpopRst$pop1,ydata$population)]
meltpopRst$region2<-ydata$ecoregion[match(meltpopRst$pop2,ydata$population)]
meltpopRst$group1<-ydata$perugroup[match(meltpopRst$pop1,ydata$population)]
meltpopRst$group2<-ydata$perugroup[match(meltpopRst$pop2,ydata$population)]
meltpopRst$lang1<-ydata$Language_Clustering[match(meltpopRst$pop1,ydata$population)]
meltpopRst$lang2<-ydata$Language_Clustering[match(meltpopRst$pop2,ydata$population)]



# Boxplots and bar plots -----------------------------
library(ggplot2)
library(ggpubr)
library(ggridges) 
### Pairwise distances *********************************************

## boxplots of pairwise distances
g<-ggplot(meltpop3)+   #all north coast comparisons including outsiders as whole
  geom_boxplot(aes(x = reorder(groupCouple, averageDifference, FUN = median), y = averageDifference))+
  coord_flip()+
  theme_pubr()
g
ggsave("boxAverageDiffPerugroup.pdf", useDingbats = FALSE)

g<-ggplot(filter(meltpop3, group1 == "Coast"| group2 == "Coast"))+   #only coast comparisons
  geom_boxplot(aes(x = groupCouple, y = averageDifference))+
  coord_flip()+
  theme_pubr()
g
ggsave("boxAverageDiffCoast.pdf", useDingbats = FALSE)

g<-ggplot(meltpop3, na.rm = TRUE)+   #region comparisons
  geom_boxplot(aes(x = reorder(regionCouple, averageDifference, FUN = median), y = averageDifference))+
  coord_flip()+
  theme_pubr()
g


# neighbour joining trees --------------------------------------
###function to make a neighbour joining tree

library(ape)
library(stats)
library(ggpubr)
library(RColorBrewer)
makeNJ = function(mat, grouping, popTable, pdfname, selection){ 
  # calculate tree
  mat = mat[which(rownames(mat)%in%selection),which(rownames(mat)%in%selection)]
  m6 = as.dist(mat, diag = F, upper = F)
  treeNJ<-nj(m6)
  treeNJ$edge.length[treeNJ$edge.length < 0] = 0.002 # manual correct for negative edges
  # get matching pops from table
  rownames(popTable)<-popTable$Population
  popTable = popTable[popTable$Population %in% rownames(mat),] #get correct populations for average distances tree
  popTableSort = data.frame(Population = rownames(mat))
  popSorted = plyr::join(popTableSort, popTable, by = "Population") # put by = "Population"!
  # make colours for groups
  groups<-table(popSorted[[grouping]])
  if(grouping == "MacroGroup"){
    colorcod = coloursMacroG[1:nlevels(popSorted[[grouping]])]
  }else if(grouping == "perugroup"){
    colorcod = coloursPeruG[1:nlevels(popSorted[[grouping]])]
  }else if(grouping == "Population"){
    colorcod = wes_palette("GrandBudapest1",4)
  }
  #colorcod = brewer.pal(n = nlevels(popSorted[[grouping]]), name = "Dark2") #enter desired grouping variable
  #colorcod<-rainbow(length(region))
  #colorcod=c("lightgreen","plum4","darkorange","brown","coral3","darkgreen", "lightblue","pink","blue","deeppink3", "gray50", "darkgoldenrod2", "magenta", "red", "gray10")
  colorstring<-c(rep(NA, nrow(popSorted)))
  for (i in 1:nrow(mat)) {
    colorstring[which(popSorted[[grouping]]==labels(groups[i]))]<-colorcod[i] #enter desired grouping variable
  }
  ## plot treeNJ
  pdf(pdfname ,paper='special')
  plot.phylo(treeNJ, type="u", tip.col=colorstring, cex=0.4)
  unlist(names(groups))->nomilegenda
  #legend("topleft",legend = nomilegenda,text.col=colorcod, cex = 0.6)
  dev.off()
}
makeNJ(matRst1, grouping = "MacroGroup", popTable, "NJ_Rstclose17_new.pdf", closerst) #groupings: "MacroGroup", "perugroup"


##############################
# plot the lines of distance on a map-----
### like in Barbieri 2017 Scientific Reports Chachapoyas

library(maps)
library(geosphere)
library(rworldmap)

## sharing map using TMRCA data (get it from appropriate section), this works

meltpopTMRCA = filter(meltpopTMRCA, pop1 != "Chogo" & pop2 != "Chogo") # because Chogo has no coordinates

meltpopTMRCA2<-meltpopTMRCA[which(meltpopTMRCA$identical>median(meltpopTMRCA$identical[meltpopTMRCA$identical!=0])),]# to clean up a bit i use a threshold above the median - excluding zeros

meltpopTMRCAca = filter(meltpopTMRCA2,(meltpopTMRCA2$pop1 %in% filter(popTable, perugroup != "OUT")$Population)&
                          (meltpopTMRCA2$pop2 %in% filter(popTable, perugroup != "OUT")$Population))
meltpop3<-meltpop2[which(meltpop2$TMRCA1000>median(meltpop2$TMRCA1000[meltpop2$TMRCA1000!=0])),]
meltpop3<-meltpop2[which(meltpop2$TMRCA2000>median(meltpop2$TMRCA2000[meltpop2$TMRCA2000!=0])),]

col.1 <- adjustcolor("darkorange", alpha=0.3)
col.2 <- adjustcolor("darkred", alpha=0.3)
edge.col <- colorRampPalette(c(col.1, col.2), alpha = TRUE)(100)

setwd("~/switchdrive/NorthCoastPeru_Master/Ychromosome/Rplots")
pdf("MapSharingIdentical4.pdf",paper='special')
pdf("MapSharingTMRCA500_3.pdf",paper='special')
pdf("MapSharingTMRCA1000_3.pdf",paper='special')
pdf("MapSharingTMRCA2000_3.pdf",paper='special')

map(database = "world", regions = ".",xlim=c(-90,-60),ylim=c(-20,12), col="grey90", fill=TRUE,  lwd=0.1)
points(x=filter(popTable, perugroup != "OUT" )$lon, y=filter(popTable, perugroup != "OUT" )$lat, pch=19,  cex=0.3, col="grey24")

for(i in 1:nrow(meltpopTMRCAca))  {
  node1 <- popTable[popTable$Population == as.character(meltpopTMRCAca[i,]$pop1),]
  node2 <- popTable[popTable$Population == as.character(meltpopTMRCAca[i,]$pop2),]
  
  arc <- gcIntermediate( c(node1[1,]$lon, node1[1,]$lat), 
                         c(node2[1,]$lon, node2[1,]$lat), 
                         n=1, addStartEnd=TRUE )
  edge.ind <- round(100*meltpopTMRCAca[i,]$identical / max(meltpopTMRCAca$identical))
  
  lines(arc, col=edge.col[edge.ind], lwd=edge.ind/10)
}
text(x=popTable$lon, y=popTable$lat, labels=popTable$code, pch=19,  cex=0.1, col="white")
#title(main = "Sharing of identical haplotypes")
dev.off()


# variance and haplotype diversity ------------------------
### original script from Barbieri et al. 2014, modified by Lea in May 2022
## don't need to run again
## variance per pop
spop <- table(ydata$population) # number of individuals for each population
spop<-spop[which(spop>4)]			#trick to have populations n>4)
popTable<-popTable[which(popTable[,1]%in%unlist(labels(spop))),]
npop <- length(spop) # number of populations
popN <- names(spop) # names of the populations
ydataRED<-ydata[which(ydata$population%in%popN),]


sequenzone<- vector("list",npop )
sequenzoneCH<- vector("list",npop )
for (i in 1:npop) {
  sequenzone[[i]]<-ydataRED[(which(ydataRED$population==popN[i])),9:25]
}
for (i in 1:npop) {
  sequenzoneCH[[i]]<-apply(sequenzone[[i]],1,paste, collapse="")
}
resDiversity<- matrix(0,ncol=4,nrow=npop,dimnames=list(popN, c("variance","popSize", "nHaplotypes", "HDiversity")) )
for (i in 1:npop) {
  resDiversity[i,1]<- mean(apply(sequenzone[[i]], 2,function(x){var(x, na.rm=T)}))
}
for (i in 1:npop) {
  resDiversity[i,2]<- spop[i]
}
for (i in 1:npop) {
  resDiversity[i,3]<- length(table(sequenzoneCH[[i]]))
}
for (i in 1:npop) {
  a<-length(table(sequenzoneCH[[i]]))
  b<-c()
  for (j in 1:a){
    b[j]<-(table(sequenzoneCH[[i]])[j]/spop[i])^2
  }
  resDiversity[i,4]<-(1-sum(b))*(spop[i]/(spop[i]-1))
}
write.table(resDiversity, "diversitySTR.txt", sep="\t")


resDiversityDF = tibble(pop = rownames(resDiversity), variance = resDiversity[,1], HDiversity = resDiversity[,4])

ggplot(data = subset(resDiversityDF, !is.na(variance)))+
  geom_bar(stat = "identity", aes(x = reorder(pop, variance), y = variance)) + coord_flip()+
  theme_pubr(base_size = 4)
ggsave("VarianceOrdered.pdf")

## maps to visualize diversity and variance
library(maps)
library(geosphere)
library(rworldmap)
library(ggplot2)
# put lat and lon information into resDiversityDF
popTlonlat = select(popTableArctic, Population, lon, lat, MacroGroup, code)
popTlonlat = rename(popTlonlat, pop = Population)
resDiversityDF = plyr::join(popTlonlat, resDiversityDF, by = "pop")
resDiversityDF = resDiversityDF[!is.na(resDiversityDF$HDiversity),] # some pops don't have value

# map with circle size proportional to HDiversity
# put colour of points into df
colorcod = coloursMacroG[1:nlevels(resDiversityDF$MacroGroup)]
region<-table(resDiversityDF[["MacroGroup"]])
colorstring<-c(rep(NA, length(resDiversityDF[["MacroGroup"]])))
for (i in 1:nrow(resDiversityDF)) {
  colorstring[which(resDiversityDF[["MacroGroup"]]==labels(region[i]))]<-colorcod[i]
}
resDiversityDF$colour = colorstring

## map with colour gradient depending on variance / diversity
# with ggplot2
map.world <- map_data(map="world")
gg <- ggplot()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=c(-47,70), xlim=c(-140,-45))   ## coordinates of lat and lon to frame the map area to plot whole america
gg2<-gg+coord_quickmap(ylim=c(-20,11), xlim=c(-85,-65))  #for central andes region and up to panama

# make palette
bluepalette = colorRampPalette(brewer.pal(8, "Blues"))(10)
bluepalette2 = colorRampPalette(colors = c("#56B1F7","#132B43"))(10)
gg + 
  geom_point(data=resDiversityDF, aes(x=lon, y=lat, colour=variance), size=1.5) +
  geom_point(data=filter(resDiversityDF, pop %in% coastpops), aes(x=lon, y=lat, colour=variance), size=2, shape = "square") +
  scale_colour_continuous(type = "gradient", low = "#56B1F7", high = "#132B43")+
  geom_text(data=resDiversityDF, aes(x=lon, y=lat, label=code), size=0.5, color = "white")+
  theme(text = element_text(size=12))+
  theme_pubr(legend = "right")
ggsave("MapVarianceGradientAll.pdf", useDingbats = FALSE)
# same for diversity 
gg + 
  geom_point(data=resDiversityDF, aes(x=lon, y=lat, colour=HDiversity), size=1.5) +
  geom_point(data=filter(resDiversityDF, pop %in% coastpops), aes(x=lon, y=lat, colour=HDiversity), size=2, shape = "square") +
  #scale_colour_distiller(palette = "Blues" , direction = 1 , values = c(0, 0.5, 0.9, 0.95, 1))+
  scale_colour_gradientn(colours = bluepalette2, values = c(0, 0.5, 0.9, 0.95, 0.99, 1))+
  #geom_point(data=filter(resDiversityDF, HDiversity == 1), aes(x=lon, y=lat, colour="purple"), size=1) +
  geom_text(data=resDiversityDF, aes(x=lon, y=lat, label=code), size=0.5, color = "white")+
  theme(text = element_text(size=12))+
  theme_pubr(legend = "right")
ggsave("MapDiversityGradientAll3.pdf", useDingbats = FALSE)
