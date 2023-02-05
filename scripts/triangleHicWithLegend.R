## Load libraries
library(plotgardener)
library(InteractionSet)

## Define common params
p <- pgParams(assembly = "hg38",
              
              ## Define positions
              x = 0.5,
              y = 0.5,
              width = 4,
              height = 2,
              
              ## Define region of interest
              chrom = "7",
              chromstart = 54660000,
              chromend = 55350000,
              
              ## Hi-C parameters
              resolution = 10e03,
              zrange = c(0, 200),
              norm = "SCALE")


## Hi-C plot
png(file="plots/triangleHicWithLegend.png",
    width = 4.05, height = 2.05, res=300, units="in")

pageCreate(width=5, height=3, showGuides=FALSE)

plot <- plotHicTriangle(data = "data/raw/hic/MEGA_K562_WT_0_inter.hic",
                        params = p)

annoHeatmapLegend(plot=plot,
                  params=p,
                  x = p$x + 0.1,
                  width=0.1,
                  height=p$height*0.5,
                  fontcolor="black")

dev.off()