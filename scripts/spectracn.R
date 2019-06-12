library(ggplot2)
library(scales)

setwd("/Users/rhiea/Genome10k/spectracn/bTaeGut2")

x_max=250
y_max=25000000
coverage=read.table("bTaeGut2.spectra-hap.hist", header=T)
head(coverage)

## Line graph
ggplot(data=coverage, aes(x=kmer_multiplicity, y=Count, group=kmer, colour=kmer)) +
  geom_line() +
  theme_bw() + scale_color_brewer(palette = "Set1") +
  scale_y_continuous(labels=comma) +
  coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))

## Area under the curve filled
ggplot(data=coverage, aes(x=kmer_multiplicity, y=Count)) +
  geom_ribbon(aes(ymin=0, ymax=pmax(Count,0), fill=kmer, colour=kmer), alpha=0.5, linetype=1) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", direction=1) +  # , breaks=rev(levels(coverage$kmer))
  scale_fill_brewer(palette = "Set1", direction=1) +   # , breaks=rev(levels(coverage$kmer))
  scale_y_continuous(labels=comma) +
  coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))

coverage=read.table("bTaeGut2.spectra-cn.hist", header=T)
head(coverage)
### Switch the order
levels(coverage$kmer)
coverage$kmer=factor(coverage$kmer, levels=levels(coverage$kmer)[c(6,2,3,4,5,1)])
levels(coverage$kmer)

## Spectra-cn, line graph
ggplot(data=coverage, aes(x=kmer_multiplicity, y=Count, group=kmer, colour=kmer)) +
  geom_line() +
  theme_bw() + scale_color_brewer(palette = "Set1") +
  scale_y_continuous(labels=comma) +
  coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))

## Area under the curve filled
ggplot(data=coverage, aes(x=kmer_multiplicity, y=Count)) +
  geom_ribbon(aes(ymin=0, ymax=pmax(Count,0), fill=kmer, colour=kmer), alpha=0.5, linetype=1) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", direction=1) +  # , breaks=rev(levels(coverage$kmer))
  scale_fill_brewer(palette = "Set1", direction=1) +   # , breaks=rev(levels(coverage$kmer))
  scale_y_continuous(labels=comma) +
  coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))

## Stacked, area filled
coverage$kmer=factor(coverage$kmer, levels=levels(coverage$kmer)[c(6,5,4,3,2,1)])
levels(coverage$kmer)
ggplot(data=coverage, aes(x=kmer_multiplicity, y=Count, group=kmer, fill=kmer)) +
  geom_area(size=.2, alpha=.8) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", direction=-1, breaks=rev(levels(coverage$kmer))) +
  scale_fill_brewer(palette="Set1", direction=-1, breaks=rev(levels(coverage$kmer))) +
  scale_y_continuous(labels=comma) +
  coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max)) # xlim(0, x_max) + ylim(0,y_max)
