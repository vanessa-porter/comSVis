.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))

### -------------------------------------------------------------------------------
### READ IN FILES
### -------------------------------------------------------------------------------

cn <- read.delim("/projects/hpv_nanopore_prj/htmcp/comSVis/output/scratch/HTMCP-03-06-02260-event1/region_cna.txt", header = T)
bedpe <- read.delim("/projects/hpv_nanopore_prj/htmcp/comSVis/output/scratch/HTMCP-03-06-02260-event1/subset.bedpe", header = F)

# Make help options
option_list = list(
  make_option(c("-b", "--bedpe"), type="character", default=NULL,
              help="subsetted bedpe file", metavar="character"),
  make_option(c("-c", "--cna"), type="character", default=NULL,
              help="CN of SV segments", metavar="character"),
  make_option(c("-d", "--depth"), type="character",
              help="Region depth bed file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default = "regionCN.txt",
              help="Output file name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$out

### -------------------------------------------------------------------------------
### MAKE THE COORDINATE DF
### -------------------------------------------------------------------------------
new_cn <- cn
new_cn$cn_change_category[grep('HPV',new_cn$chr)] <- ifelse(new_cn$cn_rounded[grep('HPV',new_cn$chr)] > 0, "gain", "no_change")
new_cn$cn_change[grep('HPV',new_cn$chr)] <- new_cn$cn_rounded[grep('HPV',new_cn$chr)]
new_cn$cn_change[grep('HPV',new_cn$chr)] <- ifelse(new_cn$cn_change[grep('HPV',new_cn$chr)] == 0,NA,new_cn$cn_change[grep('HPV',new_cn$chr)])

new_cn$posRemap <- NA
new_cn$endRemap <- NA

bedpe <- bedpe[,-c(2,5)]

chrs <- unique(c(bedpe$V1, bedpe$V4))

## lines
bedpe$x1 <- NA
bedpe$x2 <- NA
bedpe$y1 <- NA
bedpe$y2 <- NA

for (chr in chrs) {
  start = min(new_cn$pos[new_cn$chr == chr])
  
  # remap to a new coordinate
  new_cn$posRemap[new_cn$chr == chr] <- (new_cn$pos[new_cn$chr == chr] - start) + 1
  new_cn$endRemap[new_cn$chr == chr] <- (new_cn$end[new_cn$chr == chr] - start) + 1
  
  bedpe$x1[bedpe$V1 == chr] <- (bedpe$V3[bedpe$V1 == chr] - start) + 1
  bedpe$x2[bedpe$V4 == chr] <- (bedpe$V6[bedpe$V4 == chr] - start) + 1
}


### size of the whole length
sizes <- new_cn %>%
  group_by(chr) %>%
  summarise(length = max(end) - min(pos))
size <- sum(sizes$length)

### divide by 10 to get the spacing
space <- size / 15

for (n in 2:length(chrs)) {
  chr <- chrs[n]
  end <- max(new_cn$endRemap[new_cn$chr == chrs[n-1]])
  
  new_cn$posRemap[new_cn$chr == chr] <- (new_cn$posRemap[new_cn$chr == chr] + end) + space
  new_cn$endRemap[new_cn$chr == chr] <- (new_cn$endRemap[new_cn$chr == chr] + end) + space
  
  new_sv$posRemap[new_sv$V1 == chr] <- (new_sv$posRemap[new_sv$V1 == chr] + end) + space
  new_sv$endRemap[new_sv$V4 == chr] <- (new_sv$endRemap[new_sv$V4 == chr] + end) + space
  
  bedpe$x1[bedpe$V1 == chr] <- (bedpe$x1[bedpe$V1 == chr] + end) + space
  bedpe$x2[bedpe$V4 == chr] <- (bedpe$x2[bedpe$V4 == chr] + end) + space
}

new_cn$cn_rounded <- ifelse(new_cn$cn_rounded == 0, NA, new_cn$cn_rounded)


### -------------------------------------------------------------------------------
### BOXES
### -------------------------------------------------------------------------------

max = max(new_cn$cn_rounded[complete.cases(new_cn$cn_rounded)])

topedge = max / 50

sideedge = size/50

box <- data.frame(chr = chrs, ymin =  rep(0-topedge, length(chrs)), ymax = rep(topedge+max, length(chrs)), xmin = NA, xmax = NA)

for (i in 1:length(chrs)) {
  start <- min(new_cn$posRemap[new_cn$chr == chrs[i]])
  end <- max(new_cn$endRemap[new_cn$chr == chrs[i]])
  
  box$xmin[i] <- start - sideedge
  box$xmax[i] <- end + sideedge
}

### -------------------------------------------------------------------------------
### SV CONNECTIONS
### -------------------------------------------------------------------------------

## lines
new_cn <- new_cn[new_cn$pos != new_cn$end,]

for (i in 1:nrow(bedpe)) {
  p1 <- bedpe$V3[i]
  p2 <- bedpe$V6[i] 
  strand1 <- bedpe$V9[i] 
  strand2 <- bedpe$V10[i] 
  if (strand1 == "-" & strand2 == "+"){
    bedpe$y1[i] <- new_cn$cn_rounded[new_cn$pos == p1]
    bedpe$y2[i] <- new_cn$cn_rounded[new_cn$end == p2] 
  } else if (strand1 == "+" & strand2 == "-"){
    bedpe$y1[i] <- new_cn$cn_rounded[new_cn$end == p1]
    bedpe$y2[i] <- new_cn$cn_rounded[new_cn$pos == p2] 
  } else if (strand1 == "+" & strand2 == "+"){
    bedpe$y1[i] <- new_cn$cn_rounded[new_cn$end == p1]
    bedpe$y2[i] <- new_cn$cn_rounded[new_cn$end == p2] 
  } else {
    bedpe$y1[i] <- new_cn$cn_rounded[new_cn$pos == p1]
    bedpe$y2[i] <- new_cn$cn_rounded[new_cn$pos == p2] 
  } 
}

# translocations

tra <- bedpe[bedpe$V11 == "TRA",]
dup <- bedpe[bedpe$V11 == "DUP",]
del <- bedpe[bedpe$V11 == "DEL",]
inv <- bedpe[bedpe$V11 == "INV",]


### -------------------------------------------------------------------------------
### PLOT
### -------------------------------------------------------------------------------

plot <- ggplot()+ 
  # box around chromosomes
  geom_rect(data=box, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#ededed")  +
  # outer regions
  geom_rect(data=new_cn, 
          aes(xmin = posRemap, xmax = endRemap, ymin = (cn_rounded - 0.3), ymax = (cn_rounded + 0.3), fill = cn_change_category, colour = cn_change_category), size = 0.5)  +
  geom_segment(data = tra, aes(y = y1, yend = y2, x = x1, xend = x2), linetype = 2, size = 0.5)+
  geom_curve(data = dup, aes(y = y1, yend = y2, x = x1, xend = x2), linetype = 2, size = 0.5, curvature = -0.5, color = "#a0025c") +
  geom_curve(data = del, aes(y = y1, yend = y2, x = x1, xend = x2), linetype = 2, size = 0.5, curvature = 0.5) +
  geom_curve(data = inv, aes(y = y1, yend = y2, x = x1, xend = x2), linetype = 2, size = 0.5, curvature = 0.5, color = "#6D67E4") +
  theme_void() + 
  labs(x = "# of basepairs", y = "copy number") +
  scale_fill_brewer(palette = "Set1") + 
  scale_colour_brewer(palette = "Set1") + 
  theme(axis.text = element_text(),
        axis.title.y = element_text(face = "bold", size = 12, angle = 90),
        axis.title.x = element_text(face = "bold", size = 12))

plot


