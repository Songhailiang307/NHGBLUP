#########################################################################
# 								
# Function: Genetic effects of traits including additive effects and dominant effects were simulated. 									
# Written by Hailiang Song
# Copyright:  Beijing Academy of Agriculture and Forestry Sciences. All rights reserved.
# email: songhl0317@163.com									
# First version: 2022-12-6 
# 								
#########################################################################

#' @examples
#' ## Example: ploidy=4 and dominance degree=1;
#' Rscript simFInaAD.R 4 0.3 1

library(AlphaSimR)
library(dplyr)


Args <- commandArgs(trailingOnly=TRUE)
ploidy=as.numeric(Args[1]) # ploidy num
h_2=as.numeric(Args[2])     # heritability
v_d=as.numeric(Args[3])     # dominance variance


# Creating Founder population
nInd      = 600
nChr      = 10
segSites  = 1000
nploidy = ploidy

# founder population haplotypes
#set.seed(29121983)
founderPop = runMacs2(nInd,
                     nChr,
                     segSites,
                     inbred = FALSE,
                     ploidy = nploidy,
                     Ne = 50,
                     bp = 2e+07,
                     genLen = 1,
                     mutRate = 2e-6)

SP   = SimParam$ 
  new(founderPop)$                           # Initialize SimParam
  addTraitAD(250,meanDD=v_d,varDD=0.2)$        # Set QTL number
  addSnpChip(750)$                           # Set SNP number (+ QTL =total number)
  setVarE(h2=h_2)$                           # Set heritability
  setSexes("yes_sys")                        # Set gender usage


SP$quadProb = 0.15                           # the probability of quadrivalent pairing in an autopolyploid


# Modeling the Breeding Program
# Initialize the population
pop       = newPop(founderPop,simParam=SP)        #  Create initial population
population=NULL
poph2=popVD=NULL
for(cycle in 1:10){                   #  20 generations of breeding
  # Select and mate
  print(cycle)
  pop = selectCross(pop = pop, 
                    nFemale = 100,
                    nMale = 50,
                    use="rand",       #  random mate
                    nProgeny=20,      #  number of offspring
                    nCrosses=100,     #  number of crosses
                    balance = TRUE,
                    simParam=SP)
  popdata<-data.frame(pop@id,cycle,pop@father,pop@mother,pop@sex,pop@pheno,pop@gv)
  population<-rbind(population,popdata)
  poph2[cycle]=varA(pop)/varP(pop)
  popVD[cycle]=varD(pop)
  if(cycle==9){
    G9<-pullSegSiteGeno(pop)
  }else if (cycle==10){
    G10<-pullSegSiteGeno(pop)
  }
}

colnames(population)<-c("ID","G","Sire","Dam","Sex","Phe","Tbv")
print(mean(poph2))
print(mean(popVD))

write.table(population,"population.txt",quote = F,row.names =F,col.names =T)
write.table(G9,"G9.txt",quote = F,row.names =F,col.names =F)
write.table(G10,"G10.txt",quote = F,row.names =F,col.names =F)
