## Load libraries
library(plotgardener)
library(InteractionSet)

## Define common params
p <- pgParams(assembly = "hg38",
              
              ## Define region of interest
              chrom = "7",
              chromstart = 54660000,
              chromend = 55350000,
              
              ## Hi-C parameters
              resolution = 5e03,
              zrange = c(0, 60),
              norm = "SCALE")

## Hi-C plot
png(file="plots/rectangleHic.png", width = 4, height = 3, res=300, units="in")
plotHicRectangle(data = "data/raw/hic/MEGA_K562_WT_0_inter.hic",
                 params = p)
dev.off()