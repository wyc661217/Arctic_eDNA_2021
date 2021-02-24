
##########################
#This R script for subset the bam file headers to output a new header section containing headers that only appeard in alignments
#Yucheng Wang wyc661217@gmail.com

#You need first to generate 2 header text files with samtools:
#samtools view $bam_file -H > header_full.txt
#samtools view  $bam_file | cut -f3 | sort -u > header_new.txt

#then run this script with the 2 header text files as input


args <- commandArgs(TRUE)
header_full <- args[1]
header_new <- args[2]

df1 = read.csv(header_new, header = F, stringsAsFactors = F)
df2 = unique(df1[,1])

header = read.csv(header_full, header = F, stringsAsFactors = F, sep = "\t")
header = header[which(header[,1] %in% df2),]
header$V1 = paste("SN:",header$V1,sep = "")
header$v3 = "@SQ"
header = header[,c(3,1,2)]

write.table(header, "header_subset.txt", row.names = F, col.names = F, quote = F,sep = "\t")

