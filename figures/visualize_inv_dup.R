library(devtools)
library(SVbyEye)
library(stringr)
library(GenomicRanges)
#remove.packages("SVbyEye")
library(devtools)
#devtools::install_github("daewoooo/SVbyEye",force = TRUE)
library(SVbyEye)
library(stringr)
library(GenomicRanges)
library(randomcoloR)

setwd("Marmoset_Genes/Alignment")




Pafs <- Sys.glob(file.path("*/Query.progressive.100b.paf"))

for (Paf in Pafs){
  
  Dir <- unlist(str_split(Paf,"/",2))[1]
  Chrom = Dir
  paf.table <- readPaf(paf.file = Paf, include.paf.tags = TRUE, restrict.paf.tags = 'cg')
  
  Order <- unique(c(paf.table$t.name,paf.table$q.name))

  
  GeneDf <- read.table(paste0(Dir,"/Gene.bed"),sep='\t',header=F)[,c(1:3)]
  colnames(GeneDf) <- c('seqnames', 'start', 'end')

  GeneDf.gr <- makeGRangesFromDataFrame(GeneDf, keep.extra.columns = TRUE)
  GeneDf.gr$Gene <- Dir
  GeneDf.gr$y.id <- as.character(seqnames(GeneDf.gr))

  plt <- plotAVA(paf.table = paf.table,
                 seqnames.order=Order,
                 color.by = 'direction',
                 color.palette = c('+' = 'royalblue', '-'= 'gold'),
                 outline.alignments = F)
  
  plt <- addAnnotation(ggplot.obj = plt, annot.gr = GeneDf.gr, shape = 'rectangle', y.label.id = 'y.id',
                       annotation.level = 0.0, fill.by = 'Gene', coordinate.space = "self",
                       color.palette = c(Dir = 'black'))

  
  png(paste0(Dir,"_aln.png"), width = 4000, height = 2500, res=300)
  print(plt)
  dev.off()
  
  pdf(paste0(Dir,"_aln.pdf"), width = 15, height = 12)
  print(plt)
  dev.off()

}
