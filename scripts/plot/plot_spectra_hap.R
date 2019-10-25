#!/usr/bin/env Rscript

require("argparse")
require("ggplot2")
require("scales")

parser <- ArgumentParser(description = "Make spectra-hap plots. Line and filled pectra-hap plots will be generated.")
parser$add_argument("-f", "--file", type="character", help=".spectra-hap.hist file (required)", default=NULL)
parser$add_argument("-o", "--output", type="character", help="output prefix (required)")
parser$add_argument("-x", "--xdim", type="double", default=6, help="width of plot [default %(default)s]")
parser$add_argument("-y", "--ydim", type="double", default=3, help="height of plot [default %(default)s]")
parser$add_argument("-m", "--max", type="integer", default=150, help="maximum limit for k-mer multiplicity [default %(default)s]")
args <- parser$parse_args()

spectra_hap_plot  <-  function(hist, name, w=6, h=3, x_max=150) {
  
  dat=read.table(hist, header = TRUE)
  y_max=max(dat[dat[,1]!="read-total" & dat[,1]!="read-only",]$Count)
  y_max=y_max*1.5
  print(paste("y_max:", y_max, sep=" "))

  dat[,1]=factor(dat[,1], levels=unique(dat[,1]), ordered=TRUE) # Lock in the order  
  
  ## Line graph
  ggplot(data=dat, aes(x=kmer_multiplicity, y=Count, group=dat[,1], colour=dat[,1])) +
    geom_line() +
    theme_bw() +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(labels=comma) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max)) +
    labs(fill="k-mer", colour="k-mer")
  ggsave(file = paste(name, 'hap.ln.png', sep = "."), height = h, width = w)
  ggsave(file = paste(name, 'hap.ln.pdf', sep = "."), height = h, width = w)
  
  ## Area under the curve filled
  ggplot(data=dat, aes(x=kmer_multiplicity, y=Count)) +
    geom_ribbon(aes(ymin=0, ymax=pmax(Count,0), fill=dat[,1], colour=dat[,1]), alpha=0.5, linetype=1) +
    theme_bw() +
    scale_color_brewer(palette = "Set1", direction=1) +  # , breaks=rev(levels(dat$Assembly))
    scale_fill_brewer(palette = "Set1", direction=1) +   # , breaks=rev(levels(dat$Assembly))
    scale_y_continuous(labels=comma) +
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max)) +
    labs(fill="k-mer", colour="k-mer")
  ggsave(file = paste(name, 'hap.fl.png', sep = "."), height = h, width = w)
  ggsave(file = paste(name, 'hap.fl.pdf', sep = "."), height = h, width = w)
}

spectra_hap_plot(args$file, args$output, args$xdim, args$ydim, args$max)
