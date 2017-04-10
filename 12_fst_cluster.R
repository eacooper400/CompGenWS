### Week 12: Scaling Up; Part Two
### Fst Re-visited

setwd("/home/ecoope4/CompGenWS/")
source("LC_functions.R")

args=commandArgs(trailingOnly = TRUE)
input=args[1]
output=args[2]

sampleVCF <- my.read.vcf(file=input, comment.char="", header=TRUE, stringsAsFactors=FALSE)

f.column=grep("FORMAT", colnames(sampleVCF))
gen.start=f.column + 1
gen.end=length(colnames(sampleVCF))
pop1.index = grep("^FU_", colnames(sampleVCF)[gen.start:gen.end])
pop2.index = grep("^MA_", colnames(sampleVCF)[gen.start:gen.end])

fst.results <- data.frame()

for (i in 1:nrow(sampleVCF)) {
    vcf.row <- as.vector(sampleVCF[i,], mode="character")
    my.genotypes <- get.genotypes(vcf.row=vcf.row, start=gen.start, end=gen.end, format=vcf.row[f.column])
    ht <- get.Hexp(my.genotypes)
    hexp1 <- get.Hexp(my.genotypes[pop1.index])
    hexp2 <- get.Hexp(my.genotypes[pop2.index])
    N <- sum(count.genotypes(my.genotypes)) 
    n1 <- sum(count.genotypes(my.genotypes[pop1.index]))  
    n2 <- length(grep("0|1", my.genotypes[pop2.index])) 
    hs <- ((n1/N) * hexp1) + ((n2/N) * hexp2)
    fst <- (ht - hs)/ht
    fst.results = rbind(fst.results, c(vcf.row[1], vcf.row[2], ht, hexp1, hexp2, fst))
}

colnames(fst.results) <- c("CHR", "POS", "Het_T", "Het_FU", "Het_MA", "FST")

### Create a chromosome plot to look at the spatial distribution of Fst along the chromosome
pdf(file=output, width=8.5, height=11)
plot(fst.results$POS, fst.results$FST, xlab="Chromosomal Position (bp)", ylab=expression(F[ST]), pch=20, main="Spatial Distribution of FST", col="hotpink3")

### Add a line showing where mean Fst is:
abline(h=mean(fst.results$FST))

### Add another line to show where the 99th percentile is
q <- quantile(fst.results$FST, probs=c(0.9,0.95,0.99))
abline(h=q["99%"], lty=2)

### Finally, add a leged to the plot
legend(1.25e07, 0.8, c("Mean", "99th Percentile"), lty=c(1,2), cex=0.7)
dev.off()

  
