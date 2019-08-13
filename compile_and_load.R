source("C:/Users/beve/Documents/gllvm-1/R/gllvm.auxiliary.R")
source("C:/Users/beve/Documents/gllvm-1/R/gllvm.TMB.R")

library(TMB)
compile(file="C:/Users/beve/Documents/gllvm-1/src/gllvm2.cpp")
dyn.unload(dynlib("C:/Users/beve/Documents/gllvm-1/src/gllvm2"))
dyn.load(dynlib("C:/Users/beve/Documents/gllvm-1/src/gllvm2"))
library(mvabund)
data(spider)
y<-spider$abund
y[y>0]<-1
gllvm.TMB(y,num.lv=2,family="binomial",row.eff="random") #diagonal structure doesn't work atm but C code is correct.

#two
compile(file="C:/Users/beve/Documents/gllvm-1/src/gllvm2.cpp")
dyn.unload(dynlib("C:/Users/beve/Documents/gllvm-1/src/gllvm2"))
dyn.load(dynlib("C:/Users/beve/Documents/gllvm-1/src/gllvm2"))
gllvm.TMB(y,family="binomial")

#for ordinal
library(vegan)
data(dune)
y<-data.matrix(dune)
set.seed(2)
y<-matrix(sample(1:3,30*30,replace=T),ncol=30,nrow=30)

gllvm.TMB(y,num.lv=2,family="ordinal") #diagonal structure doesn't work atm but C code is correct.

##notes:
#find max value in the matrix
#iterate through levels and populate an array so I dont need to do that on the R side
#look at francis' code how to do it in C