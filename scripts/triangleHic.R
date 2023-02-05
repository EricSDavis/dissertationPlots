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
              resolution = 10e03,
              zrange = c(0, 200),
              norm = "SCALE")

## Hi-C plot
png(file="plots/triangleHic.png", width = 4, height = 2, res=300, units="in")
plotHicTriangle(data = "data/raw/hic/MEGA_K562_WT_0_inter.hic",
                params = p)
dev.off()