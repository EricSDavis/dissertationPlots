## Load libraries
library(plotgardener)
library(InteractionSet)

initialHic <- function(p = p) {
  ## Hi-C plots
  hic <- 
    plotHicSquare(data = "data/raw/hic/MEGA_K562_WT_4320_inter.hic",
                  params = p,
                  half = "both",
                  draw = FALSE)
  return(hic)
}

## Define loop
loop <- 
  GInteractions(GRanges("7:54750000-54825000"),
                GRanges("7:55225000-55300000"))

## Define function for plotting region
plotRegion <- function(p = p, hic = hic) {
  
  ## Hi-C plots
  hic2 <- 
    plotHicSquare(data = "data/raw/hic/MEGA_K562_WT_4320_inter.hic",
                  params = p,
                  half = "both")
  
  ## Annotate heatmap legend
  annoHeatmapLegend(plot = hic2,
                    params = p,
                    orientation = "v",
                    fontcolor = "black",
                    x = p$x + p$width + p$space,
                    width = p$space,
                    height = p$height*0.5)
  
  ## Label genome
  annoGenomeLabel(plot = hic2,
                  params = p,
                  y = p$y + p$height,
                  scale = "Mb")
  
  ## Label genome
  annoGenomeLabel(plot = hic2,
                  params = p,
                  axis = "y",
                  x = p$x,
                  just = c("right", "top"),
                  scale = "Mb")
  
}

addLoop <- function(hic=hic, loop=loop) {
  
  ## Annotate loop
  annoPixels(plot = hic,
             data = loop,
             shift = 0)
  
}

addLoopHighlight <- function(p=p, hic=hic, loop=loop) {
  
  annoHighlight(plot = hic,
                params = p,
                chrom = as.character(seqnames(anchors(loop, 'first'))),
                chromstart = start(anchors(loop, 'first')),
                chromend = end(anchors(loop, 'first')))
  annoHighlight(plot = hic,
                params = p,
                chrom = as.character(seqnames(anchors(loop, 'second'))),
                chromstart = start(anchors(loop, 'second')),
                chromend = end(anchors(loop, 'second')))
  
}

addOffDiagonal <- function(p=p, hic=hic, loop=loop) {
  plotHicSquare(data = "data/raw/hic/MEGA_K562_WT_4320_inter.hic",
                params = p,
                chrom = as.character(seqnames(anchors(loop, 'first'))),
                chromstart = start(anchors(loop, 'first')),
                chromend = end(anchors(loop, 'first')),
                altchrom = as.character(seqnames(anchors(loop, 'second'))),
                altchromstart = start(anchors(loop, 'second')),
                altchromend = end(anchors(loop, 'second')),
                half = "bottom")
}


## Initialize page
pageCreate(width = 8, height = 4, showGuides = FALSE)

## Define common params
p <- pgParams(assembly = "hg38",
              
              ## Define region of interest
              chrom = "7",
              chromstart = 54660000,
              chromend = 55350000,
              
              ## Hi-C parameters
              resolution = 25e03,
              zrange = c(0, 1400),
              norm = "SCALE",
              
              ## Define reference point
              x = 0.5,
              y = 0.5,
              width = 3,
              height = 3,
              length = 3,
              space = 0.1)

## Initial Hi-C object
hic <- initialHic(p = p)

## Plot region
png(file = "plots/hicMapExample01.png", width = 8, height = 4, res=300, units='in')
pageCreate(width = 8, height = 4, showGuides = FALSE)
plotRegion(p = p, hic = hic)
dev.off()

## Show genomic loci and loop
png(file = "plots/hicMapExample02.png", width = 8, height = 4, res=300, units='in')
pageCreate(width = 8, height = 4, showGuides = FALSE)
plotRegion(p = p, hic = hic)
addLoop(hic = hic, loop = loop)
# addLoopHighlight(p = p, hic = hic, loop = loop)
dev.off()

## Show off diagonal zoom of loop
p2 <- p
p2$x <- 4.0
p2$y <- 1.0
p2$width <- 1.0
p2$height <- 1.0
png(file = "plots/hicMapExample03.png", width = 8, height = 4, res=300, units='in')
pageCreate(width = 8, height = 4, showGuides = FALSE)
plotRegion(p = p, hic = hic)
addLoop(hic = hic, loop = loop)
addOffDiagonal(p = p2, hic = hic, loop = loop)
plotText(label = "25-Kb Pixels",
         x = p2$x,
         y = p2$y - p2$space/4,
         just = c('left', "bottom"))
dev.off()

## Finer resolution (10Kb)
png(file = "plots/hicMapExample04.png", width = 8, height = 4, res=300, units='in')
pageCreate(width = 8, height = 4, showGuides = FALSE)
p3 <- p
p3$resolution <- 10e03
p3$zrange <- c(0, 300)
plotRegion(p = p3, hic = hic)
addLoop(hic = hic, loop = loop)

## Show off diagonal zoom of loop
p2$x <- 4.0
p2$y <- 1.0
p2$width <- 1.0
p2$height <- 1.0
addOffDiagonal(p = p2, hic = hic, loop = loop)
plotText(label = "25-Kb Pixels",
         x = p2$x,
         y = p2$y - p2$space/4,
         just = c('left', "bottom"))

p3$x <- 5.25
p3$y <- 1.0
p3$width <- 1.0
p3$height <- 1.0
addOffDiagonal(p = p3, hic = hic, loop = loop)
plotText(label = "10-Kb Pixels",
         x = p3$x,
         y = p3$y - p3$space/4,
         just = c('left', "bottom"))
dev.off()

## Finer resolution (5Kb)
png(file = "plots/hicMapExample05.png", width = 8, height = 4, res=300, units='in')
pageCreate(width = 8, height = 4, showGuides = FALSE)
p4 <- p
p4$resolution <- 5e03
p4$zrange <- c(0, 100)
plotRegion(p = p4, hic = hic)
addLoop(hic = hic, loop = loop)

## Show off diagonal zoom of loop
p2$x <- 4.0
p2$y <- 1.0
p2$width <- 1.0
p2$height <- 1.0
addOffDiagonal(p = p2, hic = hic, loop = loop)
plotText(label = "25-Kb Pixels",
         x = p2$x,
         y = p2$y - p2$space/4,
         just = c('left', "bottom"))

p3$x <- 5.25
p3$y <- 1.0
p3$width <- 1.0
p3$height <- 1.0
addOffDiagonal(p = p3, hic = hic, loop = loop)
plotText(label = "10-Kb Pixels",
         x = p3$x,
         y = p3$y - p3$space/4,
         just = c('left', "bottom"))

p4$x <- 6.5
p4$y <- 1.0
p4$width <- 1.0
p4$height <- 1.0
addOffDiagonal(p = p4, hic = hic, loop = loop)
plotText(label = "5-Kb Pixels",
         x = p4$x,
         y = p4$y - p4$space/4,
         just = c('left', "bottom"))
dev.off()
