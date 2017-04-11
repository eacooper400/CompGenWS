### Read in the LC_functions.R file

source("/home/ecoope4/CompGenWS/LC_functions.R")

args=commandArgs(trailingOnly = TRUE)
input=args[1]
output=args[2]

sampleVCF <- my.read.vcf(file=input, comment.char="", header=TRUE, stringsAsFactors=FALSE)

out.con=file(output, open="wt")
sink(out.con, append=TRUE, type="output")

header=paste(c("Distance","D","Rsquared"), sep="\t")
cat(header, "\n")

format.col <- grep("FORMAT", colnames(sampleVCF))

for (i in 1:(nrow(sampleVCF)-1)) {
    row1=as.vector(sampleVCF[i,], mode="character") 
    g1=get.genotypes(row1, start=(format.col+1), end=ncol(sampleVCF), format=row1[format.col]) 
    pA=allele.freq(count.genotypes(g1))[1] 
    for (j in (i+1):nrow(sampleVCF)) { 
        row2=as.vector(sampleVCF[j,], mode="character")
        g2=get.genotypes(row2, start=(format.col+1), end=ncol(sampleVCF), format=row2[format.col])
        pB=allele.freq(count.genotypes(g2))[1]
        haps=as.vector(sapply(1:length(g2), FUN=get.haplotypes, genotypes1=g1, genotypes=g2))
        pAB=length(grep("00", haps))/length(haps)
        distance=as.numeric(row2[2])-as.numeric(row1[2])
        D=pAB-(pA*pB)
        rsq=(D**2)/(pA*(1-pA)*pB*(1-pB)) 
        LD.results=paste(c(distance,D,rsq), sep="\t")
        cat(LD.results, "\n")
    }
}
sink()
close(out.con)


