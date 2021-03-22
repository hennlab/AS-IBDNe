# AUTHOR: SSG
#
# This script accepts a plink .bim file and a genetic map file that the bim file should be trimmed to fit
# This script returns list of SNPs to remove from the plink file
# The required arguments are: Input .bim prefix (the .bim file to be trimmed)
#                           : Input genetic map file
#
# NOTE - this script expects the genetic map to have a header and for the first column to be the physical position
#
# USAGE:
# $ Rscript trim_bimfile.R ${FILE}_chr21 genetic_map_chr21_combined_b37.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Two arguments must be supplied", call.=FALSE)
} else if (length(args)>2) {
  stop("Too many arguments", call.=FALSE)
}

BIM <- read.table(paste(args[1], "bim", sep="."), stringsAsFactors=F)
GENMAP <- read.table(args[2], stringsAsFactors=F, header=T)
OUTPUT <- paste(args[1], "removeSNPs", sep="_")

map.range <- range(GENMAP[,1])
SNPs.to.remove <- BIM[BIM$V4 < map.range[1] | BIM$V4 > map.range[2],]

write.table(SNPs.to.remove, OUTPUT, quote=F, row.names=F, col.names=F, sep=" ")
