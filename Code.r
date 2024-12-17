library(GenVisR)
# Download file: http://genomedata.org/gen-viz-workshop/GenVisR/ALL1_CaptureDepth.tsv
# load the ouput of samtools depth into R
seqData <- read.delim("ALL1_CaptureDepth.tsv", header=F)
colnames(seqData) <- c("chr", "position", "samp1", "samp2", "samp3", "samp4")
seqData <- seqData[,c(3:6)]

# Count the occurrences of each coverage value
# install.packages("plyr")
library(plyr)
seqCovList <- apply(seqData, 2, count)

# rename the columns for each dataframe
renameCol <- function(x){
    colnames(x) <- c("coverage", "freq")
    return(x)
}
seqCovList <- lapply(seqCovList, renameCol)

# create framework data frame with entries for the min to the max coverage
maximum <- max(unlist(lapply(seqCovList, function(x) max(x$coverage))))
minimum <- min(unlist(lapply(seqCovList, function(x) min(x$coverage))))
covFramework <- data.frame("coverage"=minimum:maximum)

# Merge the framework data frame with the coverage
# seqCovList <- lapply(seqCovList, function(x, y) merge(x, y, by="coverage", all=TRUE), covFramework)
seqCovList <- lapply(seqCovList, merge, covFramework, by="coverage", all=TRUE)

# merge all data frames together
seqCovDataframe <- Reduce(function(...) merge(..., by="coverage", all=T), seqCovList)

# set all NA values to 0
seqCovDataframe[is.na(seqCovDataframe)] <- 0

# set the rownames, remove the extra column, and convert to a matrix
rownames(seqCovDataframe) <- seqCovDataframe$coverage
seqCovDataframe$coverage <- NULL
seqCovMatrix <- as.matrix(seqCovDataframe)

# rename columns
colnames(seqCovMatrix) <- c("sample1", "sample2", "sample3", "sample4")

# run covBars
covBars(seqCovMatrix)

# ceiling pileups to 1200
column_sums <- colSums(seqCovMatrix[1200:nrow(seqCovMatrix),])
column_sums <- t(as.matrix(column_sums))
rownames(column_sums) <- 1200
seqCovMatrix2 <- seqCovMatrix[1:1199,]
seqCovMatrix2 <- rbind(seqCovMatrix2, column_sums)

# run covBars
covBars(seqCovMatrix2)

# change the colours in our Plots
colorRamp <- rainbow(1200)[1:1050]
covBars(seqCovMatrix2, colour=colorRamp)

