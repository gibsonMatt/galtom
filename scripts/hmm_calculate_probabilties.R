#Matt Gibson
#2019
#Indiana University Departement of Biology
#Moyle Lab


#This file contains functions and code used to calculate emission probabilities from pariwise divergence data.
#Code for calculating these probabilities is shown only for population MG114. Identical pipelines were used for MG116, MG117, and MG115.

#setwd("/Users/matthew/Box Sync/Projects/Galapagos/Sequencing/pipeline2/introgression_analysis")

library(tidyverse)
library(readr)
library(plotly)
library(gridExtra)
library(adegenet)

#################################################
#Utility Functions
#################################################

#Binomial model 1
log_prob_pi <- function(x){
  s <- min(round(as.numeric(x[3])), round(as.numeric(x[5])))
  med115 <- round(as.numeric(x[6]))
  return(dbinom(med115, s, as.numeric(x[4]), log=T))
}

#Binomial model 2
log_prob_dxy <- function(x){
  s <- min(round(as.numeric(x[3])), round(as.numeric(x[5])))
  med115 <- round(as.numeric(x[7]))
  return(dbinom(med115, s, as.numeric(x[2]), log=T))
}

#Binomial model 3
log_prob_het <- function(x){
  s <- min(round(as.numeric(x[3])), round(as.numeric(x[5])))
  medhet <- round(as.numeric(x[8]))
  return(dbinom(medhet, s, as.numeric(x[9]), log=T))
}

#This prepares the input files for the HMM
prepHMM_focal114 <- function(chr, mg115Ind){
  #Pi and dxy to all
  dxy <- windowed.pi %>% filter(Chr == chr,a == mg115Ind, str_detect(b, '115')) %>%group_by(window) %>%summarise(X = median(pi), S = median(sites))
  pi <- windowed.pi %>% filter(Chr == chr, a == mg115Ind, str_detect(b, '114')) %>%group_by(window)  %>%summarise(X = median(pi), S = median(sites))

  #medians
  medians_to114 <- windowed.pi %>% filter(Chr==chr, a==mg115Ind, str_detect(b, '114')) %>% group_by(window, b) %>% 
  summarise(d = diffs) %>% group_by(window)%>% summarise(med = median(d))
  medians_to115 <- windowed.pi %>% filter(Chr==chr, a==mg115Ind, str_detect(b, '115')) %>% group_by(window, b) %>% 
  summarise(d = diffs) %>% group_by(window)%>% summarise(med = median(d))
  
  
  probs <- as_tibble(inner_join(dxy, pi, by="window"))
  probs <- inner_join(probs, medians_to114, by="window")
  probs <- inner_join(probs, medians_to115, by="window")
  names(probs) <- c('window', "dxy", "dxy_sites", "pi", "pi_sites", "med114", "med115")

  probs$medHet <- (probs$med114+probs$med115)/2

  probs$het <- (probs$dxy+probs$pi)/2
  print(probs)
  probs$P_114 <- apply(probs, 1, FUN = log_prob_pi)
  probs$P_115 <- apply(probs, 1, FUN = log_prob_dxy)
  probs$P_het <- apply(probs, 1, FUN = log_prob_het)
  
  write.table(probs, paste('./', 'hmmdata/', chr, '_', mg115Ind,'_MG115_het.txt', sep=''), sep='\t')
  return(probs)
}

#Plots HMM results for debugging
plot_HMM_mg114 <- function(chr, mg115Ind){
  f <- paste('./hmmdata/', chr, '_', mg115Ind, '_MG115.txt.tsv', sep = "")
  #par(mfrow=c(2,1))
  d1 <- read.csv(f, sep ='\t')
  d1$case <- ifelse(d1$annt == '+', 'gray60', funky(5)[1])
  return(d1)
  palette(c('gray60', funky(5)[1]))
  #plot.new()
  plot(d1$dxy, type='l', xlab="Window (100kb)", ylab='dXY', col="black", main=paste("Chr", chr, ": Focal ", mg115Ind, sep =''))
  for (i in seq(1, dim(d1)[1])){
    r <- d1[i,]
    abline(v = i, col = alpha(r$case, 0.2), lwd=3)
  }
  d1$case <- ifelse(d1$annt == '+', 'gray60', funky(5)[1])
  palette(c('gray60', funky(5)[1]))
  #plot.new()
  plot(d1$pi, type='l', xlab="Window (100kb)", ylab='pi (MG115)', col="black", main="")
  for (i in seq(1, dim(d1)[1])){
    r <- d1[i,]
    abline(v = i, col = alpha(r$case, 0.2), lwd=3)
  }
}

intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

splt <- function(X, pattern, col, to.order = FALSE){
  str_split(string = X, pattern = pattern, simplify = TRUE)[,col] 
}



#################################################
#Read in and pre-process pairwise divergence data generated with stacks-pairwise.py
#################################################
put.admx <-  c("MG115.3", "MG115.5", "MG115.7", "MG115.16", "MG115.17R", "MG115.20") #this is a list of individuals at site MG115 which are known to be hybrids
x <- read.csv("./mg115_mg114_pairwise_pi.txt",sep = "\t",as.is = TRUE, header = TRUE) %>%
  filter(!Chr %in% c("SL3.0ch00", "")) %>%
  dplyr::select(-LocusID)                     %>%
  mutate(Chr = str_remove(string = Chr, pattern = "SL3.0ch"))

uniques <- which(!duplicated(sapply( str_split(names(x),"_"), function(X){paste(sort(X), collapse = "_")})))

#some tidyverse magic
x <- x%>% 
  dplyr::select(uniques) %>%
  gather(key = comp, value = obs, -Chr, -StartPos) %>%
  filter(  !(is.na(obs) | obs == "") )%>%
  mutate(diffs  = splt(X = obs,  pattern = "/", col = 1),
         sites  = splt(X = obs,  pattern = "/", col = 2), 
         a      = splt(X = comp, pattern = "_", col = 1, to.order = TRUE),
         b      = splt(X = comp, pattern = "_", col = 2))        %>%
  dplyr::select(-obs, - comp)                                           %>%
  mutate_at(.vars = c("StartPos", "diffs", "sites"), as.numeric) %>%
  filter(a != b)                                                 %>% 
  mutate_at(.vars = c("a","b"), function(X){ifelse(X == "MG115.19","MG114.xx",X)})  %>%
  mutate(num114   = as.numeric(str_detect(a,"114")) + as.numeric(str_detect(b,"114")),
         num.put  = as.numeric(a %in% put.admx) + as.numeric(b %in% put.admx),
         type     = case_when(num114 == 2 ~ "114",
                              num114 == 0 & num.put == 0 ~ "115",
                              num114 == 1 & num.put == 1 ~ "put x 114",
                              num114 == 0 & num.put == 1 ~ "put x 115",
                              num.put == 2 ~ "put x put",
                              num114 == 1 & num.put == 0 ~ "114 x 115" ))


#################################################
#Convert divergence estimates per-radtag to 100kb window estimates
#################################################
window.size <- 1e5
windowed.pi <- x  %>% 
  mutate(window = cut(StartPos,breaks = seq(from = 0, 
                                            to = max(StartPos) + window.size,
                                            by =  window.size)))%>%
  group_by(a,b,type, window, Chr) %>% 
  summarise(diffs = sum(diffs),
            sites = sum(sites),
            pi = sum(diffs) / sum(sites)) %>%
  ungroup() %>% 
  mutate(type = fct_reorder(type,pi,mean)) %>% 
  mutate(pi1 = cut(pi,breaks = c(-.01,0,0.001, 0.0025 ,0.005,0.01, 0.02,0.04,1)))

windowed.pi <- windowed.pi %>% filter(!a %in% put.admx, !b %in% put.admx)
chrs <- as.character(levels(as.factor(windowed.pi$Chr)))



#################################################
#Calculate emission probabilities for each individual and write to files
#################################################

#HMM input files are written to a temp directory and are then input into introgression_hmm.py
for (c in chrs){
    prepHMM_focal114(c, 'MG114.9')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.12')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.8')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.4')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.11R')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.5')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.7')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.10')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.xx')
}

for (c in chrs){
   prepHMM_focal114(c, 'MG114.15')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.6')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.1')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.13')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.14')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.2')
}

for (c in chrs){
    prepHMM_focal114(c, 'MG114.3')
}



