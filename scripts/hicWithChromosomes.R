## Load libraries
library(plotgardener)
library(grid)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## Define chromosomes (where A < B)
chrA <- 6
chrB <- 7

## Fractional height
chrRatio <- seqlengths(seqinfo(txdb))[chrB]/seqlengths(seqinfo(txdb))[chrA]

## Initialize page
png(file = "plots/hicWithChromosomes.png", width = 8, height=4, res=300, units='in')
pageCreate(width = 6.3, height = 4, showGuides = FALSE)

## Define common params
p <- pgParams(assembly = "hg38",
              
              ## Define region of interest
              chrom = as.character(chrA),
              chromstart = 1,
              chromend = seqlengths(seqinfo(txdb))[chrA],
              
              ## Hi-C parameters
              resolution = 250e03,
              zrange = c(0, 600),
              norm = "NONE",
              
              ## Define reference point
              x = 0.5,
              y = 0.5,
              width = 3,
              height = 3,
              length = 3,
              space = 0.1)


## Chromosome A
hicA <- 
  plotHicSquare(data = "data/raw/hic/MEGA_K562_WT_4320_inter.hic",
                params = p,
                half = "both",
                just = c("left", "top"))

## Chromosomes A and B
inter <-
  plotHicSquare(data = "data/raw/hic/MEGA_K562_WT_4320_inter.hic",
                params = p,
                chrom = as.character(chrA),
                chromend = seqlengths(seqinfo(txdb))[chrA],
                altchrom = as.character(chrB),
                altchromstart = 1,
                altchromend = seqlengths(seqinfo(txdb))[chrB],
                zrange = c(0, 100),
                x = p$x + p$width,
                y = 0.5,
                half = "bottom",
                just = c("left", "top"),
                draw = FALSE)

## Customize width (since chromosome B is shorter)
inter$grobs$vp$width <- unit(p$width*chrRatio, "inches")

name <- inter$grobs$childrenOrder[2]

inter$grobs$children[[name]]$y <-
  inter$grobs$children[[name]]$y*chrRatio

inter$grobs$children[[name]]$x <-
  inter$grobs$children[[name]]$x/chrRatio

## Draw chromosome A x B
grid.draw(inter$grobs)

## Ideograms
plotIdeogram(chrom = paste0("chr", chrA),
             assembly = "hg38",
             orientation = 'h',
             x = p$x,
             y = p$y + p$height + p$space/2,
             width = p$width,
             height = p$space*2)

## Flip verical ideogram 180 degrees
vIdeogram <- 
  plotIdeogram(chrom = paste0("chr", chrA),
               assembly = "hg38",
               orientation = 'v',
               x = p$x - p$space/2 - p$space*2,
               y = p$y + p$height,
               width = p$space*2,
               height = p$height,
               just = c("right", "top"),
               draw = FALSE)
vIdeogram$grob$vp$angle <- 90
grid.draw(vIdeogram$grob)

plotIdeogram(chrom = paste0("chr", chrB),
             assembly = "hg38",
             orientation = 'h',
             x = p$x + p$width,
             y = p$y + p$height + p$space/2,
             width = p$width*chrRatio,
             height = p$space*2)

## Annotate text
plotText(label = paste("Chromosome", chrA),
         x = p$x,
         y = p$y + p$height + p$space/2 + p$space*3,
         just = c("left", "top"),
         fontcolor = "grey30")

plotText(label = paste("Chromosome", chrA),
         x = p$x - p$space/2 - p$space*3,
         y = p$y + p$height,
         just = c("left", "bottom"),
         rot = 90,
         fontcolor = "grey30")

plotText(label = paste("Chromosome", chrB),
         x = p$x + p$width,
         y = p$y + p$height + p$space/2 + p$space*3,
         just = c("left", "top"),
         fontcolor = "grey30")

## Annotate legends
annoHeatmapLegend(plot = hicA,
                  orientation = 'h',
                  fontcolor = "grey30",
                  x = p$x,
                  y = p$y - p$space,
                  width = p$width*0.75,
                  height = p$space/2)

annoHeatmapLegend(plot = inter,
                  orientation = 'h',
                  fontcolor = "grey30",
                  x = p$x + p$width,
                  y = p$y - p$space,
                  width = p$width*chrRatio*0.75,
                  height = p$space/2)
dev.off()
