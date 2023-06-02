#!/usr/bin/env Rscript
# install.packages('ggseqlogo', repos="https://mirrors.nju.edu.cn/R/")
# install.packages('seqinr', repos="https://mirrors.nju.edu.cn/R/")
# install.packages("ggplot2", repos="https://mirrors.nju.edu.cn/R/")




library(ggplot2)
library(ggseqlogo)
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)
name=args[1]
a1 <- read.fasta(paste0('./result/',name,'.fa'))
temp_str = vector(mode = 'character')
for (j in 1:length(a1)){
    temp_str[j] = toupper(c2s(a1[[j]]))
}
p <- ggseqlogo(temp_str) +
    ggtitle(name)+
    theme(
        axis.text.x = element_text(size = 4,face = "bold")
    )
ggsave(paste0('./result/',name,'.png'),plot=p, width = 13.3, height = 3.725 ,dpi = "retina")