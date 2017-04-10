### Read a VCF file into a data frame (Week 1)
my.read.vcf <- function(file, special.char="##", ...) {
    my.search.term <- paste0(special.char, ".*")  # Making a search term that looks like: "##.*", tells R to find anything containing the pattern "##" followed by anything (* is wildcard)
    clean.lines <- sub(my.search.term, "", readLines(file)) # Replace any line containing the search term with nothing (in other words remove it)
    clean.lines2 <- sub("#CHROM", "CHROM", clean.lines) # Replace the #CHROM term in the header with CHROM, so R doesn't treat it as a special character
    read.table(..., text=paste(clean.lines2, collapse="\n")) # Pass the cleaned up lines to read.table
}

### Function 2: Extract the Individual fields from the Sample Genotype string (Week 2)
get.field <- function(x, name, format) { #take as arguments the genotype string for a sample, the name of the field to extract, and the value of the format column (e.g. GT:GQ)
    fields <- unlist(strsplit(x, split=":")) #splits the string into a list of separate elements wherever it encounters the ":" character; unlist coerces this into a vector
    names(fields) = unlist(strsplit(format, split=":")) # This will name each element of the above vector according to the information in the Format string
    return(fields[name]) # Finds the desired element using the name of the field
}

### Function 3: Take a row of a vcf file, and return a vector of JUST the sample genotypes 
### e.g. c("0/0", "0/1","1/1") (Week 2)
get.genotypes <- function(vcf.row, start, end, format) {   #take as arguments the vcf file row, the starting column of the genotypes, the last column of the genotypes, and any additional arguments 
    individuals <- vcf.row[start:end] #subset the vcf row to get just the sample data
    genotypes <- vapply(individuals, FUN=get.field, FUN.VALUE=character(1), USE.NAMES=FALSE, name="GT", format=format) # use vapply to run the "get.field" function on EVERY element in the vector created above
    return(genotypes)
}

### Function 4: Take a vector of individual genotypes (such as that returned by get.genotypes)
### and return a named vector of the COUNTS for each genotype  (Week 3)
count.genotypes <- function(genotypes) { # genotypes here is an object returned by "get.genotypes"
    counts <- c(length(grep("0[/|\\|]0", genotypes)), (length(grep("0[/|\\|]1", genotypes)) + length(grep("1[/|\\|]0",genotypes))), length(grep("1[/|\\|]1", genotypes))) #use grep to get all instances of each pattern, with length to get the number of instances.  My grep pattern looks more complicated because it is searching for "0/0" OR "0|0" patterns-this will be important for working with phased data later.
    names(counts) <- c("AA","Aa","aa") #assign a name to each count, to make it easy to look up later
    return(counts)
}

### Function 5: Take a named vector of genotype counts, and return a vector with 2 allele frequencies (p,q) (Week 3)
allele.freq <- function(genotype.counts) { #genotype.counts is the object produced by "count.genotypes"
    n <- genotype.counts["AA"] + genotype.counts["Aa"] + genotype.counts["aa"] # get the total sample size by summing the counts
    p <- ((2*genotype.counts["AA"]) + genotype.counts["Aa"])/(2*n) # calculate p using HWE equivalency
    q <- 1-p # Calculate q as 1-p
    freqs <- unname(c(p,q))
    return(freqs)
}

### Function 6: Calculate expected Heterozygosity for a vector of genotypes (Week 4)
get.Hexp <- function(genotypes) {
    counts <- count.genotypes(genotypes)
    p <- allele.freq(counts)[1]
    h <- 2 * p * (1-p)
    return(h)
}

### Function 7: Take 2 vectors of individual genotypes (such as that returned by get.genotypes)
### and return a vector of 2 haplotypes for a specified individual (x)
get.haplotypes <- function(x, genotypes1, genotypes2) {
    a=unlist(strsplit(genotypes1[x], split="\\|"))
    b=unlist(strsplit(genotypes2[x], split="\\|"))
    h1=paste(c(a[1], b[1]), collapse="")
    h2=paste(c(a[2], b[2]), collapse="")
    return(c(h1,h2))
}

### Function 8: Calculate Waterson's Theta for a data frame in VCF format
### Return a single value for theta.w for the whole table (not a per bp estimate)
waterson.theta <- function(vcf) {
    Sn=nrow(vcf)
    f=grep("FORMAT", colnames(vcf))
    n=(ncol(vcf))-(f+1)
    nc=2*n
    a.n=sum(1/(seq(from=1, to=(nc-1), by=1)))
    theta.w=Sn/a.n
    return(theta.w)
}

### Function 9: Calculate Pi for a data frame in VCF format
### Use the Begun (2007) formula to deal with missing data
begun.pi <- function(vcf) {
    f=grep("FORMAT", colnames(vcf))
    pi=0 
    for (i in 1:nrow(vcf)) { 
        g=get.genotypes(as.vector(vcf[i,], mode="character"), start=(f+1), end=ncol(vcf), format=vcf[i,f])
        if (length(grep("\\.", g)) > 0) { g=g[-(grep("\\.", g))] }
        alleles=unlist(strsplit(g, split="[/||]"))
        j=min(c(length(grep("1", alleles)), length(grep("0", alleles))))
        n=length(alleles)
        x=((2*j)*(n-j))/(n*(n-1)) 
        pi=pi+x
    }
    return(pi)
}

### Function 10: Calculate variance of d (for Tajima's D)
### Take the number of chromosomes and the number of SNPs as arguments
variance.d <- function(n,S) {
    a1=sum(1/(seq(from=1, to=(n-1), by=1)))
    a2=sum(1/((seq(from=1, to=(n-1), by=1))**2))
    b1=(n+1)/(3*(n-1))
    b2=(2*((n**2)+n+3))/((9*n)*(n-1))
    c1=b1 - (1/a1)
    c2=b2-((n+2)/(a1*n)) + (a2/(a1**2))
    e1=c1/a1
    e2=c2/((a1**2)+a2)
    var=(e1*S) + (e2*S*(S-1))
    return(var)
}

### Function 11: Calculate Tajima's D for a data frame in VCF format
### Return a single value for D
tajimas.d <- function(vcf) {
    f=grep("FORMAT", colnames(vcf))
    k=begun.pi(vcf)
    theta=waterson.theta(vcf)
    Sn=nrow(vcf)
    n.chr=(ncol(vcf) - (f+1))*2
    var=variance.d(n=n.chr, S=Sn)
    D=(k-theta)/(sqrt(var))
    return(D)
}

### Function 12: Calculate Dxy between 2 populations
### Return a single value for Dxy
dxy=function(vcf, pop1.search, pop2.search) {
    f=grep("FORMAT", colnames(vcf))
    pop1.index=grep(pop1.search, colnames(vcf))
    pop2.index=grep(pop2.search, colnames(vcf))
    dxy=0
    for (i in 1:nrow(vcf)) { 
        g=get.genotypes(as.vector(vcf[i,], mode="character"), start=(f+1), end=ncol(vcf), format=vcf[i,f])
        if (length(grep("\\.", g)) > 0) { g=g[-(grep("\\.", g))] }
        g1=g[pop1.index]
        g2=g[pop2.index]
        alleles1=unlist(strsplit(g1, split="[/||]"))
        alleles2=unlist(strsplit(g2, split="[/||]"))
        ## Let j be the minor allele count in population 1
        j=length(grep("1", alleles1))
        ## Let k be the minor allele count in population 2
        k=length(grep("1", alleles2))
        ## Let n1 be the number of alleles in pop1, n2 the number in pop2, and N the total number of alleles
        n1=length(alleles1)
        n2=length(alleles2)
        N=n1+n2
        ## Let p1 be the probability of seeing the minor allele in population 1 AND the major allele in pop. 2
        p1 = (j/n1) * ((n2-k)/n2)
        ## Let p2 be the probability of seeing the major allele in pop. 1 AND the minor allele in pop. 2
        p2 = ((n1-j)/n1) * (k/n2)
        ## The probability that pop1 and pop2 are different at this site is the sum of p1 and p2
        x=p1 + p2
        dxy=dxy+x
    }
    return(dxy)
}
