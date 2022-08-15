#####in the postgibbs sampling file was generated running postgibbs with 0 iter and 5 thin
#column 1 is the number of saved samples usually is ordered
##column 2 is the round from where the sample was saved 
####numer of parameters  that in ther polynomias of second order 

rm(list=ls())
###checking postgibbs convergence
install.packages("boa")
library("boa") ###
setwd("C:/Users/descamil/OneDrive - purdue.edu/Purdue Folder/Thesis-Papers/Chapter 1- Multivariate analaysis/Final_documents")
data<- read.table("postgibbs_samples",header=FALSE)
#data<- data[,c(4:7)]
data<- data[,-c(1:3)]
data<- as.matrix(data)
#dimnames(link)<-list(a$V2,seq(1,132) ## when i put dinnames it didn;t work without it it runs
                     ##after this i did
#boa.menu() ## this comand oppens an interactive selection window where you select the type of analysis 
############# after checking for convergence I will extract the variance covariance matrix 

a<-colMeans(data)              

##spliting the vector a in G and R
G<- a[1:66]
library(patr1ckm)
#G<- as.matrix(G)
G<-vec2sym(G, diagonal = NULL, lower = FALSE, byrow = TRUE)
G<- as.matrix(G)
G_cor<- cov2cor(G)
dim(G_cor)
colnames(G_cor)<- c("Suc", "Raf",	"Sta", "SW", "Flo",	"Mat","Oil",	"Prot",	"Hgt",	"Yld",	"Ldg")
rownames(G_cor)<-c("Suc", "Raf",	"Sta", "SW", "Flo",	"Mat","Oil",	"Prot",	"Hgt",	"Yld",	"Ldg")


##calculating the p-value manually
t<- G_cor/sqrt((1-(G_cor)^2)/(1096-2))
t<- as.matrix(t)
p.value.matrix<- pt(abs(t), df= 1096-2, lower.tail = FALSE)
p.value.matrix<- matrix(p.value.matrix,11)
dimnames(p.value.matrix)<- list(colnames(G_cor), rownames(G_cor))

#Removing columns and rows 
G<- G[1:8,1:8]
G_cor<- G_cor[1:8,1:8]
t<- t[1:8,1:8]
p.value.matrix<- p.value.matrix[1:8,1:8]

#storing results as a maatrix


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png("Genetic_Corr_gibbs.png", width = 700, height = 600)
corrplot::corrplot(G_cor, type = "upper",method= "number", tl.srt = 50, tl.col = "black",p.mat = p.value.matrix, sig.level = 0.05, order="hclust",insig = "blank", col=cbbPalette[c(2,3)], title = "Genetic Correlations",  mar=c(0,0,4,0) )
dev.off()

##create a nice table

# Define notions for significance levels; spacing is important.
mystars <- ifelse(p.value.matrix < .0001, "****", ifelse(p.value.matrix < .001, "*** ", ifelse(p.value.matrix < .01, "**  ", ifelse(p.value.matrix < .05, "*   ", "    "))))
## trunctuate the correlation matrix to two decimal
R <- format(round(cbind(rep(-1.11, 8), G_cor), 2))[,-1]

## build a new matrix that includes the correlations with their apropriate stars
Rnew <- matrix(paste(R, mystars, sep=""), ncol=8)
diag(Rnew) <- paste(diag(R), " ", sep="")
rownames(Rnew) <- colnames(G_cor)
colnames(Rnew) <- paste(colnames(G_cor), "", sep="")
Rnew <- as.matrix(Rnew)
Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
Rnew <- as.data.frame(Rnew)

## remove last column and return the correlation matrix
Rnew <- cbind(Rnew[1:length(Rnew)-1])

write.csv(Rnew, file="gencor_gibbs.csv")


################residual var cov matrix

R<- a[67:132]
R<-vec2sym(R, diagonal = NULL, lower = FALSE, byrow = TRUE)
R_cor<- cov2cor(R)
colnames(R_cor)<- c("Suc", "Raf",	"Sta", "SW", "Flo",	"Mat",	"Oil",	"Prot",	"Hgt",	"Yld",	"Ldg")
rownames(R_cor)<-c("Suc", "Raf",	"Sta", "SW", "Flo",	"Mat",	"Oil",	"Prot",	"Hgt",	"Yld",	"Ldg")


t<- R_cor/sqrt((1-(R_cor)^2)/(1096-2))
t<- as.matrix(t)
p.value.matrix<- pt(abs(t), df= 1096-2, lower.tail = FALSE)
p.value.matrix<- matrix(p.value.matrix,11)
dimnames(p.value.matrix)<- list(colnames(R_cor), rownames(R_cor))

#Removing columns and rows 
R_cor<- R_cor[1:8,1:8]
t<- t[1:8,1:8]
p.value.matrix<- p.value.matrix[1:8,1:8]
R<-R[1:8,1:8]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png("nonGenetic_Correlationplot_gibbs.png", width = 700, height = 600)
corrplot::corrplot(R_cor, type = "upper",method= "number", tl.srt = 50, tl.col = "black",p.mat = p.value.matrix, sig.level = 0.05, order="hclust",insig = "blank", col=cbbPalette[c(2,3)], title = "Non Genetic Correlations",  mar=c(0,0,4,0) )
dev.off()

##create a nice table

# Define notions for significance levels; spacing is important.
mystars <- ifelse(p.value.matrix < .0001, "****", ifelse(p.value.matrix < .001, "*** ", ifelse(p.value.matrix < .01, "**  ", ifelse(p.value.matrix < .05, "*   ", "    "))))
## trunctuate the correlation matrix to two decimal
R <- format(round(cbind(rep(-1.11, 8), R_cor), 2))[,-1]

## build a new matrix that includes the correlations with their apropriate stars
Rnew <- matrix(paste(R, mystars, sep=""), ncol=8)
diag(Rnew) <- paste(diag(R), " ", sep="")
rownames(Rnew) <- colnames(R_cor)
colnames(Rnew) <- paste(colnames(R_cor), "", sep="")
Rnew <- as.matrix(Rnew)
Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
Rnew <- as.data.frame(Rnew)

write.csv(Rnew, file="nongencor_gibbs.csv")

##calculating heritabilites
add<- diag(G)
#add<- add[,-c(9,10,11)]
res<-diag(R)
#res<- res[,-c(9,10,11)]
names<-c("Suc", "Raf",	"Sta", "SW", "Flo",	"Mat",	"Oil",	"Prot")
table<- cbind(names, add, res)
table<- as.data.frame(table)
table$add<-as.numeric( as.character(table$add))
table$res<-as.numeric( as.character(table$res))
table<- table%>% mutate(total= add + res, h2= add/(add+res))

write.csv(table, file="narrowh2_gibbs.csv")

# library
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyr)
library(ggeasy)

table1<- table[,c(-5)]
colnames(table1)<-c("trait","Add","res", "total")
table1<- table1%>% mutate(Additive= (add)/total , Residual= (res)/total)
table1<- table1[,c(1,5,6)]
data_long <- gather(table1, component, measurement, Additive:Residual, factor_key=TRUE)
png("variancecomponents.png", width = 700, height = 600)
# Small multiple
ggplot(data_long, aes(fill=component, y=measurement, x=trait)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired")+
  ggtitle("Proportion of variance explained") +
  #theme_ipsum() +
  xlab("Traits")+
  ylab("Variance explained")+
  ggeasy::easy_center_title()+
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent"), text = element_text(size=16, face="bold"))

dev.off()

library(dplyr)
##calculating the principal components using the princomp 
colnames(G_cor)<- colnames(R_cor)
rownames(G_cor)<- rownames(R_cor)
Gen_PCA<- princomp(G_cor, cor = T)

R_PCA<- princomp(R_cor, cor = T)

##################creating PCA plot from genetic correlations

library(tidyverse)

library(factoextra)
png("PCA_addcor.png", width = 700, height = 600)

fviz_pca_var(Gen_PCA,
             repel = TRUE, col.circle = "white", col.var = "#E94B3CFF", labelsize=6) +labs(title = "Biplot of Carbohydrate and Seed Traits- Additive correlation", x= "PC1 50.7%", y= "PC2 19.2%") +theme(text = element_text(size = 17),
                                                                                                                                                                     axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()
# ##################creating PCA plot from non genetic correlations


png("PCA_Rcor.png", width = 700, height = 600)

fviz_pca_var(R_PCA,
             repel = TRUE, col.circle = "white", col.var = "#E94B3CFF", labelsize=6) +labs(title = "Biplot of Carbohydrate and Seed Traits- non-additive correlation", x= "PC1 23.8%", y= "PC2 22%") +theme(text = element_text(size = 17), axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()



#########loading phenotype######################################################################
Pheno<- read.table("Pheno1.txt")
rownames(Pheno)<- Pheno$V1
Pheno<- Pheno[,c(2:7,9:10)]

colnames(Pheno)<-  c("Suc", "Raf",	"Sta", "SW", "Flo",	"Mat",	"Oil",	"Prot")
##how to obtain the correlatin matrix and the p-values
library(dplyr)
library(Hmisc)
#Pheno1<- Pheno[,-1]
###running "Pearson" for centralized data
Scaled<-as.data.frame(scale(Pheno, center = T, scale = T))### to calculate spearman correlation do the same but changen the type
cor_1 <-rcorr(as.matrix(Scaled), type = "pearson") ### to calculate spearman correlation do the same but changen the type
p_value_pearson<- cor_1$P ## extracting p-values
corr_pearson<- cor_1$r ##extracting correlation coefficients


##calculating the principal components using the princomp 
pearson_PCA<- princomp(corr_pearson, cor = T)

##################creating PCA plot from phenotypic correlations

library(tidyverse)

library(factoextra)
png("PCA_pearsoncor.png", width = 700, height = 600)

fviz_pca_var(pearson_PCA,
             repel = TRUE, col.circle = "white", col.var = "#E94B3CFF", labelsize=6) +labs(title = "Biplot of Carbohydrate and Seed Traits- Pearson correlation", x= "PC1 42.4%", y= "PC2 19.1%") +theme(text = element_text(size = 17),
                                                                                                                                                                                                          axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()



###running "spearman"
cor_2 <-rcorr(as.matrix(Scaled), type = "spearman") ### to calculate spearman correlation do the same but changen the type
p_value_spearman<- cor_2$P ## extracting p-values
corr_Spearman<- cor_2$r ##extracting correlation coefficients

##calculating the principal components using the princomp 
spearman_PCA<- princomp(corr_Spearman, cor = T)

##################creating PCA plot from spearman correlations#################################

library(tidyverse)

library(factoextra)
png("PCA_spearmancor.png", width = 700, height = 600)

fviz_pca_var(spearman_PCA,
             repel = TRUE, col.circle = "white", col.var = "#E94B3CFF", labelsize=6) +labs(title = "Biplot of Carbohydrate and Seed Traits- Spearman correlation", x= "PC1 44.9%", y= "PC2 19.6%") +theme(text = element_text(size = 17),
                                                                                                                                                                                                         axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()
 
###creating a plot with the four correlations 

library(factoextra)
library(ggpubr)


png("PCA_plotsblack.png", width = 1200, height = 1200)

a<-fviz_pca_var(pearson_PCA,
             repel = TRUE, col.circle = "white", col.var = "black", labelsize=7) +labs(title="",x= "PC1 44.4%", y= "PC2 18.7%") +theme(text = element_text(size = 18),
                                                                                                                                                                                                         axis.title = element_text(size = 20), axis.text = element_text(size = 18))


b<-fviz_pca_var(spearman_PCA,
             repel = TRUE, col.circle = "white", col.var = "black", labelsize=7) +labs(title="",x= "PC1 48.09%", y= "PC2 18%") +theme(text = element_text(size = 18),axis.title = element_text(size = 20), axis.text = element_text(size = 18))


c<-fviz_pca_var(Gen_PCA,
             repel = TRUE, col.circle = "white", col.var = "black", labelsize=7) +labs(title="",x= "PC1 54.17%", y= "PC2 20.5%") +theme(text = element_text(size = 18),
                                                                                                                                                                                                          axis.title = element_text(size = 20), axis.text = element_text(size = 18))

d<-fviz_pca_var(R_PCA,
             repel = TRUE, col.circle = "white", col.var = "black", labelsize=7) +labs(title="", x= "PC1 27.7%", y= "PC2 25.8%") +theme(text = element_text(size = 18), axis.title = element_text(size = 20), axis.text = element_text(size = 18))


ggarrange(a, b, c,d + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2) %>% ggexport(filename= "PCA_plotsgold.png", width=800, height=800)                                                                                                                                                                                               


dev.off()



library(igraph)
ADJ_gen=huge::huge(G_cor, method='glasso',verbose=F)
plot(ADJ_gen)
ADJ_gen=huge::huge(G_cor, lambda= 0.337, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_gen)<-list(rownames(G_cor), colnames(G_cor))

png("UDGM_gencor_gibbs.png", width = 700, height = 600)
plot(igraph::graph.adjacency(adjmatrix=ADJ_gen),vertex.label.cex=2.5, arrow.mode=3, edge.curve=0.4,vertex.color="white"
     , vertex.size=14,vertex.frame.color= "black",  vertex.lable.dist=0, edge.color="black",edge.arrow.size= 0.4, frame=TRUE,
     vertex.shape= "none" )
title("Undirected Graphical Model for Genetic Correlations",cex.main=2,col.main="black")
dev.off()

ADJ_nongen=huge::huge(R_cor, method='glasso',verbose=F)
plot(ADJ_nongen)
ADJ_nongen=huge::huge(R_cor, lambda= 0.191, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_nongen)<-list(rownames(R_cor), colnames(R_cor))
png("UDGM_rescor_gibbs.png", width = 700, height = 600)
plot(igraph::graph.adjacency(adjmatrix=ADJ_nongen),vertex.label.cex=2.5, arrow.mode=3, edge.curve=0.4,vertex.color="white"
     , vertex.size=14,vertex.frame.color= "black",  vertex.lable.dist=0, edge.color="black",edge.arrow.size= 0.4, frame=TRUE,
     vertex.shape= "none" )
title("Undirected Graphical Model for Residual Correlations",cex.main=2,col.main="black")
dev.off()


###phenotypic UDGM
library(igraph)
ADJ_spearman=huge::huge(corr_Spearman, method='glasso',verbose=F)
plot(ADJ_spearman)
ADJ_spearman=huge::huge(corr_Spearman, lambda= 0.343, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_spearman)<-list(rownames(corr_Spearman), colnames(corr_Spearman))

png("UDGM_spearmancor_gibbs.png", width = 700, height = 600)
plot(igraph::graph.adjacency(adjmatrix=ADJ_spearman),vertex.label.cex=2.5, arrow.mode=3, edge.curve=0.4,vertex.color="white"
     , vertex.size=14,vertex.frame.color= "black",  vertex.lable.dist=0, edge.color="black",edge.arrow.size= 0.4, frame=TRUE,
     vertex.shape= "none" )
title("Undirected Graphical Model for Spearman Correlations",cex.main=2,col.main="black")
dev.off()

ADJ_pearson=huge::huge(corr_pearson, method='glasso',verbose=F)
plot(ADJ_pearson)
ADJ_pearson=huge::huge(corr_pearson, lambda= 0.258, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_pearson)<-list(rownames(corr_pearson), colnames(corr_pearson))
png("UDGM_pearsoncor_gibbs.png", width = 700, height = 600)
plot(igraph::graph.adjacency(adjmatrix=ADJ_pearson),vertex.label.cex=2.5, arrow.mode=3, edge.curve=0.4,vertex.color="white"
     , vertex.size=14,vertex.frame.color= "black",  vertex.lable.dist=0, edge.color="black",edge.arrow.size= 0.4, frame=TRUE,
     vertex.shape= "none" )
title("Undirected Graphical Model for Pearson Correlations",cex.main=2,col.main="black")
dev.off()




####a plot of all
ADJ_spearman=huge::huge(corr_Spearman, lambda= 0.343, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_spearman)<-list(rownames(corr_Spearman), colnames(corr_Spearman))

ADJ_pearson=huge::huge(corr_pearson, lambda= 0.301, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_pearson)<-list(rownames(corr_pearson), colnames(corr_pearson))

ADJ_gen=huge::huge(G_cor, lambda= 0.346, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_gen)<-list(rownames(G_cor), colnames(G_cor))
ADJ_nongen=huge::huge(R_cor, lambda= 0.191, method='glasso',verbose=F)$path[[1]]
dimnames(ADJ_nongen)<-list(rownames(R_cor), colnames(R_cor))

jpeg("UDGM_plots_gibbs.jpg", width=800, height=800)
par(mfrow=c(2,2),mar=c(2,2,2,2))

plot(igraph::graph.adjacency(adjmatrix=ADJ_pearson, mode = "undirected"), vertex.color="lightskyblue1"
     , vertex.size=22,vertex.frame.color= "black", arrow.mode=0, vertex.lable.dist=2, edge.color="black",vertex.label.cex=1.8,  frame=TRUE )+ title("A",cex.main=2,col.main="black")

plot(igraph::graph.adjacency(adjmatrix=ADJ_spearman, mode = "undirected"), vertex.color="lightskyblue1"
     , vertex.size=22,vertex.frame.color= "black", arrow.mode=0, vertex.lable.dist=2, edge.color="black", vertex.label.cex=1.8
    , frame=TRUE)+ title("B",cex.main=2,col.main="black")


plot(igraph::graph.adjacency(adjmatrix=ADJ_gen, mode = "undirected"), arrow.mode=0, vertex.color="lightskyblue1"
     , vertex.size=22,vertex.frame.color= "black",  vertex.lable.dist=2, edge.color="black",vertex.label.cex=1.8, 
     frame=TRUE) +title("C",cex.main=2,col.main="black")

plot(igraph::graph.adjacency(adjmatrix=ADJ_nongen, mode = "undirected"), arrow.mode=0, vertex.color="lightskyblue1"
     , vertex.size=22,vertex.frame.color= "black",  vertex.lable.dist=2, edge.color="black",vertex.label.cex=1.8, 
     frame=TRUE) +title("D",cex.main=2,col.main="black")


dev.off()
