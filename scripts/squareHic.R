## Load libraries
library(plotgardener)
library(InteractionSet)

## Define common params
p <- pgParams(assembly = "hg38",
              
              ## Define region of interest
              chrom = "7",
              chromstart = 54660000, # - 150e03,
              chromend = 55350000, # + 150e03,
              
              ## Hi-C parameters
              resolution = 10e03,
              zrange = c(0, 200),
              norm = "SCALE")

## Hi-C plot
png(file="plots/squareHic.png", width = 5, height = 5, units = 'in', res = 300)
plotHicSquare(data = "data/raw/hic/MEGA_K562_WT_0_inter.hic",
              params = p,
              half="both")
dev.off()