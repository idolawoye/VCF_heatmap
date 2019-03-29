#! /usr/bin/env Rscript

library(vcfR)
library(ape)
library(argparse)
library(RColorBrewer)

parser <- ArgumentParser(description='VCF visualisation script')
parser$add_argument("vcf", type='character', action="store", help="Input VCF file")
parser$add_argument("ref", action="store", help="Input reference fasta genome file")
parser$add_argument("gff", action="store", help="Input annotation gff file")
parser$print_help()

args <- commandArgs(trailingOnly = TRUE)

if (length(args)<3) {
  stop("Three Arguments must be supplied: VCF file, fasta file and GFF file", call.=FALSE)
} else if (length(args)==3) {
  print("Checking input files")
}

vcf_file <- args[[1]]
dna_file <- args[[2]]
annotation_file <- args[[3]]

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(annotation_file, sep="\t", quote="")

chrom <- create.chromR(name='VCF_heatmap', vcf = vcf, seq=dna, ann=gff)
genotype <- extract.gt(chrom, element = "GT", as.numeric = TRUE)
rownames(genotype) <- 1:nrow(genotype)

coul <- colorRampPalette(brewer.pal(8, "Spectral"))(25)
# Remove all rows with infinite values
genotype[!rowSums(!is.finite(genotype)),]
# Replace all non-infinite values with 0
genotype[!is.finite(genotype)] <- 0
# Create heatmap PNG output
png(file = "heatmap.png", width=960, height=960, units="px", pointsize=24)
heatmap(genotype, scale="row",col=coul, xlab = "Genomes", ylab = "SNPs")
dev.off()

noquote("DONE.")
