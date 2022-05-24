# BINARY ANALYSIS #

# read in the files
namemap=readRDS("rubysnamemap.rds")
trees=readRDS("zoonomiatrees.rds")

# make hibphenvals3 using my name map (length 173) & with vscodes2 (includes the correct names for the master tree)
hibphenvals3=namemap$pheno
names(hibphenvals3)=namemap$vscodes2
hibphenvals3=hibphenvals3[!is.na(names(hibphenvals3))] #length 173
saveRDS(hibphenvals3, "hibphenvals3.rds")

# get RERs 
hibRERS=getAllResiduals(trees, useSpecies = names(hibphenvals3))
saveRDS(hibRERS, "hibRERs.rds")

# perform the correlation with binary phenotype 
allspecs=names(hibphenvals3)
fgspecs=names(hibphenvals3)[hibphenvals3]
hibpath=foreground2Paths(fgspecs, trees, plotTree=T, clade="all", useSpecies = allspecs)
hibcors=correlateWithBinaryPhenotype(hibRERS, hibpath)
saveRDS(hibcors, "hibcors.rds")

# BAYES FACTORS #

# read in namemap
namemap = readRDS("rubysnamemap.rds")

#read in trees and RERs
trees=readRDS("zoonomiatrees.rds")
hibRER=readRDS("hibRERs.rds")

#create phenotype paths 
hibphens = readRDS("hibphenvals3.rds")
hibNoBatphens = readRDS("hibphenNoBats.rds") # for fgspecs of hibnobatpath
allspecs=names(hibphens)
# hibpath = all hibernators (including bats)
fgspecs=names(hibphens)[hibphens] 
hibpath=foreground2Paths(fgspecs, trees, plotTree=T, clade="all", useSpecies = allspecs)
# batpath = all bats (including hibernators)
fgspecs = namemap$vscodes2[namemap$isbat]
batpath=foreground2Paths(fgspecs, trees, plotTree=T, clade="all", useSpecies = allspecs)
# batnohibpath = non hibernating bats
fgspecs = namemap$vscodes2[namemap$isbat & !namemap$pheno]
batnohibpath=foreground2Paths(fgspecs, trees, plotTree=T, clade="all", useSpecies = allspecs)
# hibnobatpath = non-bat hibernators 
fgspecs = names(hibNoBatphens)[hibNoBatphens]
hibnobatpath=foreground2Paths(fgspecs, trees, plotTree=T, clade="all", useSpecies = allspecs)

#make dataframe to store Bayes factors
library(BayesFactor)
allBayesFactors=data.frame(matrix(nrow=nrow(hibRER), ncol=4))
colnames(allBayesFactors)=c("hibernator", "bat", "batnohibernator", "hibernatornobat")
rownames(allBayesFactors)=rownames(hibRER)

#calculate Bayes factors per gene
count=1
while(count<=nrow(hibRER)){
  #hibernator
  tempdf=data.frame(trait=hibpath, RER=hibRER[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$hibernator[count]=bf
  
  #bat
  tempdf=data.frame(trait=batpath, RER=hibRER[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$bat[count]=bf 
  
  #batnohibernator
  tempdf=data.frame(trait=batnohibpath, RER=hibRER[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$batnohibernator[count]=bf 
  
  #hibernatornobat
  tempdf=data.frame(trait=hibnobatpath, RER=hibRER[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$hibernatornobat[count]=bf
  
  count=count+1
  print(count)
}

saveRDS(allBayesFactors, "allBayesFactors.rds")



#permulations:

library(RERconverge)

#permulations for AA convergence GLS analyses:
tree=read.tree(file="Permulations/242m_Tree_176sps_Hib_noTorpnoMU.nh")
phen=read.table("Permulations/comparisons_242m_Phenotypes_176sps_Hib_noTorpnoMU_corrected", sep="\t", header = T)
sum(tree$tip.label %in% phen$species)
length(tree$tip.label)
phenvec=phen$pheno
phenvec=as.numeric(phenvec)
names(phenvec)=phen$species

#adjusted perm function
################################################################################
getSimBinPheno=function(mt, root, phenvec, fgnum=NULL, internal=0, drop=NULL){
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=root.phylo(mt, root, resolve.root = T)
  t=drop.tip(t, drop)
  rm=ratematrix(t, phenvec)
  sims=sim.char(t, rm, nsim = 1)
  nam=rownames(sims)
  s=as.data.frame(sims)
  simulatedvec=s[,1]
  names(simulatedvec)=nam
  top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
  # plot(t)
  return(top)
}
################################################################################

set.seed(1000)
df=data.frame(matrix(nrow=22, ncol=1000))
count=1
while(count<=1000){
  df[,count]=getSimBinPheno(tree, "Homo_sapiens", phenvec)
  count=count+1
}
saveRDS(df, "permulationfgs.rds")
write.table(df, file="Permulations/permulated_fg_vals.txt", sep="\t", quote=F, col.names = F, row.names = F)

#check species
sp=read.table("Permulations/176species.txt")
length(phenvec)
sum(tree$tip.label %in% sp$V1)
sum(sp$V1 %in% phen$species)




#RERconverge permulations:
allcors=readRDS("RubyHibernationResults/allgeneresults.rds")
bf=readRDS("RubyHibernationResults/allBayesFactors.rds")
hibcors=readRDS("RubyHibernationResults/hibcors.rds")

hibRERs=readRDS("hibRERs.rds")
trees=readRDS("zoonomiatrees.rds")
hibphenvals=readRDS("hibphenvals.rds")
hibphenvalnum=as.numeric(hibphenvals)
names(hibphenvalnum)=names(hibphenvals)

# t=foreground2Tree(names(hibphenvalnum)[hibphenvalnum==1], trees, plotTree=T, clade="all", useSpecies = names(hibphenvalnum))
# sum(t$edge.length) #35 total
# sum(hibphenvalnum) #21 tip
# sum(t$edge.length)-sum(hibphenvalnum) #14 internal

trees=trees
root="REFERENCE"
phenvec=hibphenvalnum
fgnum = 35
internal=14
drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% names(hibphenvalnum))]

statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "statdf.rds")
saveRDS(pvaldf, "pvaldf.rds")




#get perm p-val
statdf=readRDS("statdf.rds")
temp=data.frame(hibcors$Rho)
rownames(temp)=rownames(hibcors)
realandperm=merge(temp, statdf, by="row.names")
rownames(realandperm)=realandperm$Row.names
realandperm=realandperm[,-1]

permpval=data.frame(matrix(nrow=nrow(realandperm), ncol=1))
colnames(permpval)="permP"
rownames(permpval)=rownames(realandperm)

count=1
while(count<=nrow(realandperm)){
  real=realandperm[count,1]
  perms=as.numeric(realandperm[count,2:ncol(realandperm)])
  perms=perms[!is.na(perms)]
  p=sum(abs(perms)>abs(real))/length(perms)
  permpval[count, 1]=p
  print(count)
  count=count+1
}

hibcorswithpermp=hibcors
hibcorswithpermp=merge(hibcorswithpermp, permpval, by="row.names")
rownames(hibcorswithpermp)=hibcorswithpermp$Row.names
hibcorswithpermp=hibcorswithpermp[,-1]
hibcorswithpermp$permP.adj=p.adjust(hibcorswithpermp$permP, method="BH")
saveRDS(hibcorswithpermp, "hibcorswithpermp.rds")



#do more perms
allcors=readRDS("RubyHibernationResults/allgeneresults.rds")
bf=readRDS("RubyHibernationResults/allBayesFactors.rds")
hibcors=readRDS("RubyHibernationResults/hibcors.rds")

hibRERs=readRDS("hibRERs.rds")
trees=readRDS("zoonomiatrees.rds")
hibphenvals=readRDS("hibphenvals.rds")
hibphenvalnum=as.numeric(hibphenvals)
names(hibphenvalnum)=names(hibphenvals)

# t=foreground2Tree(names(hibphenvalnum)[hibphenvalnum==1], trees, plotTree=T, clade="all", useSpecies = names(hibphenvalnum))
# sum(t$edge.length) #35 total
# sum(hibphenvalnum) #21 tip
# sum(t$edge.length)-sum(hibphenvalnum) #14 internal

trees=trees
root="REFERENCE"
phenvec=hibphenvalnum
fgnum = 35
internal=14
drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% names(hibphenvalnum))]

statdf=NULL
pvaldf=NULL

count=1
while(count<=4000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "MOREstatdf.rds")
saveRDS(pvaldf, "MOREpvaldf.rds")

statdf=readRDS("statdf.rds")
pvaldf=readRDS("pvaldf.rds")
MOREstatdf=readRDS("MOREstatdf.rds")
MOREpvaldf=readRDS("MOREpvaldf.rds")
allstatdf=merge(statdf, MOREstatdf, by="row.names")
rownames(allstatdf)=allstatdf$Row.names
allstatdf=allstatdf[,-1]
allpvaldf=merge(pvaldf, MOREpvaldf, by="row.names")
rownames(allpvaldf)=allpvaldf$Row.names
allpvaldf=allpvaldf[,-1]
saveRDS(allstatdf, "allstatdf.rds")
saveRDS(allpvaldf, "allpvaldf.rds")



#get all perm pval
statdf=readRDS("allstatdf.rds")
temp=data.frame(hibcors$Rho)
rownames(temp)=rownames(hibcors)
realandperm=merge(temp, statdf, by="row.names")
rownames(realandperm)=realandperm$Row.names
realandperm=realandperm[,-1]

permpval=data.frame(matrix(nrow=nrow(realandperm), ncol=1))
colnames(permpval)="permP"
rownames(permpval)=rownames(realandperm)

count=1
while(count<=nrow(realandperm)){
  real=realandperm[count,1]
  perms=as.numeric(realandperm[count,2:ncol(realandperm)])
  perms=perms[!is.na(perms)]
  p=sum(abs(perms)>abs(real))/length(perms)
  permpval[count, 1]=p
  print(count)
  count=count+1
}

hibcorswithpermp=hibcors
hibcorswithpermp=merge(hibcorswithpermp, permpval, by="row.names")
rownames(hibcorswithpermp)=hibcorswithpermp$Row.names
hibcorswithpermp=hibcorswithpermp[,-1]
hibcorswithpermp$permP.adj=p.adjust(hibcorswithpermp$permP, method="BH")
saveRDS(hibcorswithpermp, "hibcorswithpermpALL.rds")



#run 5000 more (this should be enough)
allcors=readRDS("RubyHibernationResults/allgeneresults.rds")
bf=readRDS("RubyHibernationResults/allBayesFactors.rds")
hibcors=readRDS("RubyHibernationResults/hibcors.rds")

hibRERs=readRDS("hibRERs.rds")
trees=readRDS("zoonomiatrees.rds")
hibphenvals=readRDS("hibphenvals.rds")
hibphenvalnum=as.numeric(hibphenvals)
names(hibphenvalnum)=names(hibphenvals)

# t=foreground2Tree(names(hibphenvalnum)[hibphenvalnum==1], trees, plotTree=T, clade="all", useSpecies = names(hibphenvalnum))
# sum(t$edge.length) #35 total
# sum(hibphenvalnum) #21 tip
# sum(t$edge.length)-sum(hibphenvalnum) #14 internal

trees=trees
root="REFERENCE"
phenvec=hibphenvalnum
fgnum = 35
internal=14
drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% names(hibphenvalnum))]


#####batch 1
statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("batch 1, finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "statdf1.rds")
saveRDS(pvaldf, "pvaldf1.rds")

#####batch 2
statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("batch 2, finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "statdf2.rds")
saveRDS(pvaldf, "pvaldf2.rds")

#####batch 3
statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("batch 3, finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "statdf3.rds")
saveRDS(pvaldf, "pvaldf3.rds")

#####batch 4
statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("batch 4, finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "statdf4.rds")
saveRDS(pvaldf, "pvaldf4.rds")

#####batch 5
statdf=NULL
pvaldf=NULL

count=1
while(count<=1000){
  
  #get phenotype:
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)
  # while(blsum!=fgnum){
  while(blsum>(fgnum+5) | blsum<(fgnum-5)){
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    tf=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(tf$edge.length)
  }
  
  #get path:
  p=tree2Paths(tf, trees, useSpecies = names(hibphenvalnum))
  #run correlation:
  c=correlateWithBinaryPhenotype(hibRERs, p)
  
  if(count==1){
    statdf=data.frame(c$Rho)
    rownames(statdf)=rownames(c)
    pvaldf=data.frame(c$P)
    rownames(pvaldf)=rownames(c)
  }else{
    temp=data.frame(c$Rho)
    rownames(temp)=rownames(c)
    statdf=merge(statdf, temp, by="row.names", all = T)
    rownames(statdf)=statdf$Row.names
    statdf=statdf[,-1]
    temp=data.frame(c$P)
    rownames(temp)=rownames(c)
    pvaldf=merge(pvaldf, temp, by="row.names", all = T)
    rownames(pvaldf)=pvaldf$Row.names
    pvaldf=pvaldf[,-1]
  }
  
  print(paste0("batch 5, finished perm: ", count))
  count=count+1
}

saveRDS(statdf, "statdf5.rds")
saveRDS(pvaldf, "pvaldf5.rds")

#combine all perms

allstatdf=readRDS("allstatdf.rds")
allpvaldf=readRDS("allpvaldf.rds")

stat1=readRDS("statdf1.rds")
stat2=readRDS("statdf2.rds")
stat3=readRDS("statdf3.rds")
stat4=readRDS("statdf4.rds")
stat5=readRDS("statdf5.rds")

pval1=readRDS("pvaldf1.rds")
pval2=readRDS("pvaldf2.rds")
pval3=readRDS("pvaldf3.rds")
pval4=readRDS("pvaldf4.rds")
pval5=readRDS("pvaldf5.rds")

s=merge(allstatdf, stat1, by="row.names")
rownames(s)=s$Row.names
s=s[,-1]
s=merge(s, stat2, by="row.names")
rownames(s)=s$Row.names
s=s[,-1]
s=merge(s, stat3, by="row.names")
rownames(s)=s$Row.names
s=s[,-1]
s=merge(s, stat4, by="row.names")
rownames(s)=s$Row.names
s=s[,-1]
s=merge(s, stat5, by="row.names")
rownames(s)=s$Row.names
s=s[,-1]
saveRDS(s, "statdf10k.rds")

p=merge(allpvaldf, pval1, by="row.names")
rownames(p)=p$Row.names
p=p[,-1]
p=merge(p, pval2, by="row.names")
rownames(p)=p$Row.names
p=p[,-1]
p=merge(p, pval3, by="row.names")
rownames(p)=p$Row.names
p=p[,-1]
p=merge(p, pval4, by="row.names")
rownames(p)=p$Row.names
p=p[,-1]
p=merge(p, pval5, by="row.names")
rownames(p)=p$Row.names
p=p[,-1]
saveRDS(p, "pvaldf10k.rds")




#get all perm pval with 10k
statdf=readRDS("statdf10k.rds")
temp=data.frame(hibcors$Rho)
rownames(temp)=rownames(hibcors)
realandperm=merge(temp, statdf, by="row.names")
rownames(realandperm)=realandperm$Row.names
realandperm=realandperm[,-1]

permpval=data.frame(matrix(nrow=nrow(realandperm), ncol=1))
colnames(permpval)="permP"
rownames(permpval)=rownames(realandperm)

count=1
while(count<=nrow(realandperm)){
  real=realandperm[count,1]
  perms=as.numeric(realandperm[count,2:ncol(realandperm)])
  perms=perms[!is.na(perms)]
  p=sum(abs(perms)>abs(real))/length(perms)
  permpval[count, 1]=p
  print(count)
  count=count+1
}

hibcorswithpermp=hibcors
hibcorswithpermp=merge(hibcorswithpermp, permpval, by="row.names")
rownames(hibcorswithpermp)=hibcorswithpermp$Row.names
hibcorswithpermp=hibcorswithpermp[,-1]
hibcorswithpermp$permP.adj=p.adjust(hibcorswithpermp$permP, method="BH")
saveRDS(hibcorswithpermp, "hibcorswithpermp10k.rds")




#add BF
hibcorswithpermp=readRDS("hibcorswithpermp10k.rds")
bf=readRDS("RubyHibernationResults/allBayesFactors.rds")
justbf=bf[,1:2]
justbf$hibvsbat=justbf$hibernator/justbf$bat
colnames(justbf)=c("BayesFactorHibernator", "BayesFactorBat", "BayesFactorHibvsBat")


finalhibres=merge(hibcorswithpermp, justbf, by="row.names")
rownames(finalhibres)=finalhibres$Row.names
finalhibres=finalhibres[,-1]

saveRDS(finalhibres, "finalhibres.rds")
write.csv(finalhibres, "finalhibres.csv")


#get top genes
finalhibres=readRDS("finalhibres.rds")
ff=finalhibres[finalhibres$p.adj<0.05&finalhibres$permP.adj<0.15&finalhibres$BayesFactorHibvsBat>=5,]
ff=na.omit(ff)
saveRDS(ff, "tophibernationcors.rds")
write.csv(ff, "tophibernationcors.csv")




#Make BF plot
library(ggplot2)
library(ggpubr)

finalhibres=readRDS("finalhibres.rds")

x=finalhibres[finalhibres$p.adj<0.05,]
x$gene=rownames(x)
x$sublab=ifelse(abs(x$Rho)>0.25&x$permP.adj<0.15, x$gene, NA)
x$bfbin=ifelse(x$BayesFactorHibvsBat>5, "Hibernators", "Bats")
x$permbin=ifelse(x$permP.adj<0.15, "Significant", "Not Significant")
pdf("HibernationBFfig.pdf", width = 7, height=5)
ggplot(data=x, aes(x=Rho, y=BayesFactorHibvsBat, label=sublab, color=bfbin, shape=permbin))+
  geom_point()+
  ylim(c(0,20))+
  geom_text_repel()+
  labs(color="Bayes Factor Support", shape="Permulation p-value")+
  theme_bw()+
  ggtitle("significant genes from RERconverge")+
  xlab("Rho, non-hibernating ancestor")+
  ylab("Hibernation BF/Bat BF")
dev.off()



















