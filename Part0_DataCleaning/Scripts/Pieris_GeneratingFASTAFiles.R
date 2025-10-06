#Parsing and Characterizing Pieris Proteome

#Load Necessary Package
library(Biostrings)

#Set the correct working directory 
setwd("OneDrive\ -\ University\ of\ Kansas/PierisProteomeProject/Analysis_Results/DataCleaning/RelevantRCodeScripts/")

##Upload pierisGenome dataset
pierisGenome <- readAAStringSet(file = "../../../RawData/Pieris_rapae_proteins.fa")

#--------------------------------------------------------------------------------------------------------------
#Part 1:  Cross reference TrxID with Gene Name and Gene ID
#Update names in FASTA file using the previously generated table 
  #Output 35000 genes
##Output table

##INTERMEDIATE STEP - Used Linux: in the gff file, used the grep command in conjunction with cut to write out to file the gene names and gene IDs
#grep "\tgene\t" Pieris_rapae_pierap_codingGene_final.gff | cut -f9 > geneID.txt

##Generating the Table 
geneinfo <- read.table(file = "../../../RawData/geneID.txt", sep =";", header = FALSE, quote = "", as.is = T) #15047 gene IDs

#Update column names in the table
names(geneinfo)[1] <- "GeneID"
names(geneinfo)[2] <- "GeneName"
protNames <- names(pierisGenome)

#Remove "ID=" and "Name=" 
geneinfo$GeneID <- gsub("ID=", "", geneinfo$GeneID)
geneinfo$GeneName <- gsub("Name=", "", geneinfo$GeneName)

#Use the `match()` function to subset geneinfo$GeneID in a way that is consistent with the ordering and identity of transcript IDs on AAstringset read from fasta.
trxList <- strsplit(names(pierisGenome), split = " gene=")  #Parse fasta header from AAStringSet into list of unique IDs

prot.mtx <- do.call(what = rbind, args = trxList) #Turn this list into a matrix, to access only trxID or GeneName directly 
names(prot.mtx)[1] <- "TranscriptID"
names(prot.mtx)[2] <- "GeneName"
prot.geneID <- geneinfo$GeneID[match(prot.mtx[,2], table = geneinfo$GeneName )] #Use match to get GeneID values corresponding to GeneNames in the fasta ordering

#Combine into one table
trxID.geneName.geneID <- cbind(prot.mtx, prot.geneID)
#Fix names in table
colnames(trxID.geneName.geneID) <- c("TrxID", "GeneName", "GeneID")

write.csv(x = trxID.geneName.geneID[order(trxID.geneName.geneID[,1]), ], file = "../ResultsOutput/TrxID.GeneName.GeneID.csv", row.names = F, quote = T)

#--------------------------------------------------------------------------------------------------------------
#Part 2: Update headers in FASTA entry to contain GeneID
#Paste together to give the new IDs for the AAStringSet
## newNames <- paste(prot.mtx[,1], prot.geneID, prot.mtx[,2])
newNamesKeys <- paste0(prot.mtx[,1], " geneID=", prot.geneID, " genename=", prot.mtx[,2])  
names(pierisGenome) <- newNamesKeys
## length(names(pierisGenome)) 

#EXTRA: Clean the sequences (remove ".")
tmpchar <- as.character(pierisGenome)
tmpchar <- gsub(pattern = ".", x = tmpchar, fixed = T, replacement = "")
pierisGenome <- AAStringSet(tmpchar)

  #Write out FASTA file
writeXStringSet(x = pierisGenome, filepath = "../ResultsOutput/CleanedProteomeFASTAs/PierisRenamedProteins.fas") 

#--------------------------------------------------------------------------------------------------------------
#Part 3: Remove Redundant Sequences
#Reduce redundancy - still allows for isoforms
  #Generate Object with non-redundancy 
##Necessary Code Previously Generated 
# genesIDList <- split( x=prot.mtx[,1], f=prot.mtx[,2])
#Use split function (f = prot.mtx[,2]) 
geneAAsetsList <- split(pierisGenome, f = prot.mtx[,2]) 
#Then use lapply and unique
geneAAsetsListuniq <- lapply(geneAAsetsList, unique) #Standard list output 
## geneAAsetsListuniq <- unlist(geneAAsetsListuniq)
## str(geneAAsetsListuniq)
## length(geneAAsetsListuniq) #Sanity check

#Fix the data structure
geneAAsetsListuniq <- AAStringSetList(geneAAsetsListuniq)
names(geneAAsetsListuniq) <- NA
uniqAAstringset <- unlist(geneAAsetsListuniq) #To get back to StringSet
## names(uniqAAstringset[1:3])

#Write out FASTA file
writeXStringSet(x = uniqAAstringset, filepath = "../ResultsOutput/CleanedProteomeFASTAs/PierisNonRedundant.fas") 

#--------------------------------------------------------------------------------------------------------------
#Part 4: Reduce to longest CDS
  #Generate Object with longest CDS
#Extract the longest sequence for each gene 
#which.max gives you the index of which one is the longest or the max
#Want the length of the sequence
indexSeq <- which.max(width(geneAAsetsListuniq[[2]])) #Index that reflects the longest transcript 
## geneAAsetsListuniq[[2]][indexSeq] #Subset to get longest StringSetObject
## width(geneAAsetsListuniq[[2]][indexSeq])

#Create a function 
getLongest <- function(x) {
  indexSeq <- which.max(width(x))
  longestSSObject <- x[indexSeq]
  return(longestSSObject)
}
#In function, extract width of x and if vector seqWidths is length 1, return x
##Else run through the which.max and subset

#Subset list in lapply to only include the ones that have > 1 transcript 
## getLongest(geneAAsetsListuniq[[2]])

#Put defined function into lapply 
longestSeqList <- lapply(geneAAsetsList, FUN = getLongest) 
## longestSeqSet <- unlist(longestSeqList) 

#Fix the data structure
longestCDS <- AAStringSetList(longestSeqList) #Convert from standard list to StringSetList
## longestCDS <- unlist(longestCDS)
names(longestCDS) <- NA #Remove list names
longestCDS <- unlist(longestCDS)
##Whittle down to a single representative protein - longest one (for ortholog analysis)
#Write out FASTA file
writeXStringSet(x = longestCDS, filepath = "../ResultsOutput/CleanedProteomeFASTAs/PierisLongestCDS.fas") 

#--------------------------------------------------------------------------------------------------------------
#Part 5: Analyze the nature of redundancy

pdf(file = "../ResultsOutput/VisualizedData/GenomePlots.pdf")
#The plot below demonstrates per gene the redundant and unique versions of transcript: 
## geneID <- strsplit(names(pierisGenome), split = " gene=") 
## gene.mtx <- do.call(what = rbind, args = geneID)
## colnames(gene.mtx)[1] <- "TrxID"
## colnames(gene.mtx)[2] <- "GeneName"
## geneAAsetsList <- split(pierisGenome, f = gene.mtx[ ,2]) 
lengthSS <- sapply(geneAAsetsList, length)

## geneAAsetsListuniq <- lapply(geneAAsetsList, unique)
lengthSSUnique <- sapply(geneAAsetsListuniq, length) #Length = count of transcripts

plot(x = lengthSS, y = lengthSSUnique, pch = 19, xlab = "Original Count", 
     ylab = "Unique Transcripts", col = rgb(0,0,0,0.1))

#This barplot shows the difference (in counts) between the original and unique counts of the genome.
diffLengths <- lengthSS - lengthSSUnique
diffLengthsTable <- table(diffLengths) #Tabulate the counts

barplot(diffLengthsTable, xlab = "Count Length", ylab = "Frequency",
        main = "Difference Between Original and Unique Transcript Length", ylim = c(0,12000), space = 0)

#This barplot shows the number of unique transcripts versus gene length. 
geneLength <- sapply(longestSeqList, FUN = width) 

plot(x = geneLength, y = lengthSSUnique, pch = 19, xlab = "Gene Length", 
     ylab = "Unique Transcripts", col = rgb(0,0,0,0.1))

dev.off()
