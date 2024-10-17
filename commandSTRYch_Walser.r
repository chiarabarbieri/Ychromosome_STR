

###### Chiara Barbieri
### October 2024
### commands from the paper Resutik et al. 
# Uncovering Genetic Signatures of the Walser Migration in the Alps: Patterns of Diversity and Differentiation


#setwd("/Users/chiarabarbieri/Documents/Walser_forensics")

library(plyr)
library(tidyverse)

ydata<-read.csv("Walser_STR_YHRD_checkDYS389ii.csv", header=T, as.is=T)
ydata<-ydata[,-which(colnames(ydata)=="DYS385")]
ydata<-ydata[-which(ydata$ID=="7039")] #exclude one individual with partial haplotype

ydata$haplo<-NA
for (i in 1:nrow(ydata)){
ydata$haplo[i]<-paste(ydata[i,3:23], collapse = "_")
}

spop <- table(ydata$Location) #or ydata1
npop <- length(spop) # number of populations
popN <- names(spop) # names of the populations

## identical haplo --> verified!
countIdentical<-table(ydata$haplo, ydata$Location)
countIdentical<-as.data.frame(countIdentical)
#countIdentical$sum<-NA
sumidentical<-  apply(countIdentical, 1, sum)
countOnlyIdentical<-countIdentical[which(sumidentical>1),]
meltcoi<-melt(countOnlyIdentical)
meltcoi<-meltcoi[-which(meltcoi$value==0),]

# BRUVO distances----
library(pegas)
library(adegenet)
library(poppr)

STRgenind<-df2genind(ydata[,3:23], sep="\t", ind.names =ydata$ID, pop=ydata$Location, ploidy=1)
STRgenind$loc.n.all

#NEIdist<-nei.dist(STRgenind) # does not handle missing data and returns infinite values, takes a while

ssr <- rep(1, nLoc(STRgenind)) # assign weight


BRUVO<-bruvo.dist(STRgenind, replen = ssr) 
mBRUVO<-as.matrix(BRUVO)

write.table(mBRUVO, "mBRUVOdistmatrix_new_ssr1.txt", sep="\t")
mBRUVO = read.table("mBRUVOdistmatrix_new_ssr1.txt", sep = "\t")
mBRUVO = as.matrix(mBRUVO)
diag(mBRUVO)<-NA
#mBRUVO[upper.tri(mBRUVO, diag = TRUE)] = NA
rownames(mBRUVO)<-ydata$Location
colnames(mBRUVO)<-ydata$Location

library(reshape)
meltBruvo<-melt(mBRUVO)
meltBruvo<-meltBruvo[-which(is.na(meltBruvo$value)),]


### analysis of sharing of identical haplotypes
# IDENTIC

IDENTIC<-meltBruvo[which(meltBruvo$value==0),] 


# sharing within population
IDENTICwithin<-IDENTIC[which(IDENTIC$X1==IDENTIC$X2),]
sharingwithinpop<-table(IDENTICwithin$X1)/2 # divide by 2, the matrix was symmetrical
sharingwithinpop<-sharingwithinpop[sort(unlist(labels(sharingwithinpop)))] #sort alphabetically
FREQsharingwithinpop<-sharingwithinpop/spop

# a trick to have a double matrix but the diagonal with single values

BRUVOdiag<-mBRUVO
BRUVOdiag[upper.tri(BRUVOdiag, diag = TRUE)] = NA
meltBRUVOdiag<-melt(BRUVOdiag)
meltBRUVOdiag<-meltBRUVOdiag[-which(is.na(meltBRUVOdiag$value)),]
IDENTICdiag<-meltBRUVOdiag[which(meltBRUVOdiag$value==0),] 
IDENTICwithinremoveduplicates<-IDENTICdiag[which(IDENTICdiag$X1==IDENTICdiag$X2),]
IDENTICbetween<-IDENTIC[which(IDENTIC$X1!=IDENTIC$X2),]

IDENTICwithaTrick<-rbind(IDENTICwithinremoveduplicates,IDENTICbetween)

# sharing all population


correctlistpops<-c("Unteres  Rhonetal", "Lötschental", "Oberes Rhonetal", "Vispertal", "Törbel",
                   "Gressoney", "Bosco Gurin", "Rheinwald", "Vals", "Safiental", "Avers",
                   "Furna", "Valzeina", "Lugnez", "Thusis", "Savognin")

# tick labels: use population names together with sample size in parenthesis
correctlistpopsSIZE<-correctlistpops
  for (i in 1:length(correctlistpops)){
    correctlistpopsSIZE[i]<-paste(correctlistpops[i], " (",spop[correctlistpops][i], ")", collapse = "", sep="")
  }

spop[correctlistpops]

p<-ggplot(IDENTICwithaTrick, aes(x=X1, y=X2)) + 
  geom_count( alpha=0.8,aes( size = after_stat(n)), color="darkred") +
  guides(color = 'legend')+
  labs(x="", y="", title="") +
  #geom_point(shape=21) +
#  scale_color_gradient(low="lightblue", high="red") +  # the color code corresponds to n of haplotypes shared
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10)) +
   scale_x_discrete(limits=correctlistpops)+       #the poporder can be a different population order that you want to display in the plot
  scale_y_discrete(limits=correctlistpops, labels=correctlistpopsSIZE)
p
ggsave("Ychr_IDenticalHaplotypesSymmetricPlay_new.pdf", useDingbats=FALSE)





#******************************
#******************************


### analysis of sharing of similar haplotypes
# SIMILAR
# using thresholds of percentile distance distribution
# 0.1%, 0.5%, 1%

# repeat the analysis for each percentile threshold and save different plot names

quantt<-quantile(meltBruvo$value, probs=c(0.001))
quantt<-quantile(meltBruvo$value, probs=c(0.005))
quantt<-quantile(meltBruvo$value, probs=c(0.01))

threshold<-as.numeric(quantt)
SIMILAR<-meltBruvo[which(meltBruvo$value<=threshold),]

# a trick to have a double matrix but the diagonal with single values

BRUVOdiag<-mBRUVO
BRUVOdiag[upper.tri(BRUVOdiag, diag = TRUE)] = NA
meltBRUVOdiag<-melt(BRUVOdiag)
meltBRUVOdiag<-meltBRUVOdiag[-which(is.na(meltBRUVOdiag$value)),]
SIMILARdiag<-meltBRUVOdiag[which(meltBRUVOdiag$value<threshold),] 
SIMILARwithinremoveduplicates<-SIMILARdiag[which(SIMILARdiag$X1==SIMILARdiag$X2),]
SIMILARbetween<-SIMILAR[which(SIMILAR$X1!=SIMILAR$X2),]

SIMILARwithaTrick<-rbind(SIMILARwithinremoveduplicates,SIMILARbetween[,1:3])

# sharing all population


# correctlistpops<-c("Unteres  Rhonetal", "Lötschental", "Oberes Rhonetal", "Vispertal", "Törbel",
#                    "Gressoney", "Bosco Gurin", "Rheinwald", "Vals", "Safiental", "Avers",
#                    "Furna", "Valzeina", "Lugnez", "Thusis", "Savognin")
# 
# spop[correctlistpops]

p3<-ggplot(SIMILARwithaTrick, aes(x=X1, y=X2)) + 
  geom_count( alpha=0.8,aes( size = after_stat(n)), color="darkred") +
  guides(color = 'legend')+
  labs(x="", y="", title="") +
  #geom_point(shape=21) +
  #  scale_color_gradient(low="lightblue", high="red") +  # the color code corresponds to n of haplotypes shared
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size = 10)) +
  scale_x_discrete(limits=correctlistpops)+       #the poporder can be a different population order that you want to display in the plot
  scale_y_discrete(limits=correctlistpops, labels=correctlistpopsSIZE)
p1
p2
p3


#### ********************************
# all 4 plots together
# FIGURE S5
###

library(gridExtra)
pmulti<-grid.arrange(p, p1, p2, p3,nrow = 2)


ggsave("YchromosomeSharingProfilesNEW.pdf", pmulti, height = 8, width = 10, useDingbats=FALSE)



### distribution of BRUVO distances


library(ggridges)
theme_set(theme_minimal())

pq1<- ggplot(meltBruvo, aes(x = value, y = X1, fill = factor(after_stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE, alpha=0.7) +
  scale_fill_viridis_d(name = "Quartiles")+
  labs(x="BRUVO distances", y="", title="") 
ggsave("DistributionBruvoPops.pdf", height = 5, width = 7, useDingbats=FALSE)


meltBruvo$XX<-"_"

 pq<- ggplot(meltBruvo, aes(x = value,y=XX, fill = factor(after_stat(quantile)))) +
    stat_density_ridges(
      geom = "density_ridges_gradient", calc_ecdf = TRUE,
      quantiles = 4, quantile_lines = TRUE, alpha=0.7) +
    scale_fill_viridis_d(name = "Quartiles")+
    labs(x="BRUVO distances", y="", title="") 
  
 ggsave("DistributionBruvoTOT.pdf", height = 3, width = 7, useDingbats=FALSE)
 


 
 
 


 #******************************
 ## analysis o population relationships
 #******************************
 

# neighbour joining trees --------------------------------------
### function to make a neighbour joining tree
# mat fixfertig, grouping in quotes, name as in popTable, needs to be as factor, don't forget .pdf
library(ape)
library(stats)
library(ggpubr)
library(RColorBrewer)

# import the  matrix
# generated with YHRD tool online, RST distances
 
mat = read.table("RST_SelectionEurope2.txt", header=T, row.names = 1, sep="\t")

m6 = as.dist(mat, diag = F, upper = F)
treeNJ<-nj(m6)
treeNJ$edge.length[treeNJ$edge.length < 0] = 0.002 # manual correct for negative edges

# create the colorcode for the population names in the tree
colorstring<-rep("gray20",length(treeNJ$tip.label))
colorstring[1:15]<-"chocolate"
colorstring[c(4,14,6)]<-"blue". #romance speaking communities (excluding Lugnez)

# FIGURE S4

## plot treeNJ
pdf("NJ_selectionEurope.pdf" ,paper='special')
plot.phylo(treeNJ, type="u", tip.col=colorstring, cex=0.4)
#unlist(names(groups))->nomilegenda
#legend("topleft",legend = nomilegenda,text.col=colorcod, cex = 0.6)
dev.off()





