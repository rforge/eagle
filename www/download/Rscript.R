# Analysis of toy data sets
library(Eagle)


# read PLINK ped file
g <- ReadMarker(filename="genoDemo.ped", type="PLINK")

# read phenotypic file
p <- ReadPheno(filename="phenoDemo.dat")

# read map file
m <- ReadMap(filename="mapDemo.dat")

# perform preliminary analysis - conservative
anal <- AM(trait="trait1", fform="pc1+pc2", geno=g, pheno=p, map=m) 

# summarize results
SummaryAM(AMobj = anal, geno=g, pheno=p, map=m)



## Perform analysis with a 5% false positive rate
fp <- FPR4AM(falseposrate=0.05, trait="trait1", fform="pc1+pc2",  
             geno = g, pheno=p, map=m )

anal <- AM(trait="trait1", fform="pc1+pc2", geno=g, pheno=p, map=m, gamma=fp$setgamma) 


SummaryAM(AMobj = anal, geno=g, pheno=p, map=m)



