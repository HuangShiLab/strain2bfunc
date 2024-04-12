#!/usr/bin/env Rscript
#################################################################
# Function: UMAP
# Call: Rscript PM_Umap.R -m meta_data -d dist_file [-l T/F -o outfile -a axesfile]
# Updated at Apr. 12, 2024
# Updated by ZHANG Yufeng <yfz96@connect.hku.hk>, HUANG Shi <shihuang@hku.hk>
# Faculty of Dentistry, HKU
#################################################################

## install necessary libraries
p <- c("optparse","vegan", "ade4","ggplot2","grid","RColorBrewer","umap")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(  
  make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix file [Required]"),
  make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
  make_option(c("-l", "--drawlabel"), type="logical", default=F,help="If enable the sample label [Optional, default %default]"),
  make_option(c("-o", "--outfile"), type="character", default='Umap.pdf', help="Output UMAP coordinates [default %default]"),
  make_option(c("-p", "--pointsize"), type="double", default=7.0, help="Point size on UMAP plot [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# paramenter checking
if(is.null(opts$meta_data)) stop('Please input a meta data file')
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

axesfile<-paste(opts$outfile, ".txt", sep="")
ps <- opts$pointsize

## load data
# import meta file
meta_orig <- read.table(file=opts$meta_data, header=TRUE, row.names=1, as.is=FALSE)
meta_orig <- meta_orig[order(rownames(meta_orig)), , drop = FALSE]
# import dist file
dst_orig <- read.table(file=opts$dist_file, header=TRUE, row.names=1)
dst_orig <- dst_orig[order(rownames(dst_orig)), order(colnames(dst_orig)), drop = FALSE]

if(!identical(rownames(dst_orig), colnames(dst_orig))) {
  stop("Please ensure the consistency between the row names and column names of the distance matrix.")
}

if(!identical(rownames(dst_orig), rownames(meta_orig))) {
  stop("Please ensure the consistency between the row names of the distance matrix and that of the meta data.")
}

## main calc & draw function definition
PM_umap <- function(da, md) {
  rn <- rownames(da)
  cn <- colnames(da)
  
  meta <- md[lapply(md, class)=="factor"]
  n <- ncol(meta)
  
  sampleNumber <- length(rn)

  # dst <- as.dist(da)
  dst <- as.matrix(da)
  umap <- umap(dst, n_components = 3)
  umap <- umap$layout
  
  if(ncol(umap) < 3) {
    stop("The dimensionality of the UMAPoA plot is less than 3, making it impossible to visualize.")
  }
  
  
  colnames(umap)<-c("UMAP1","UMAP2","UMAP3")
  data<-as.data.frame(umap)
  
  # write axes file
  write.table(data,file=axesfile,sep="\t",quote=FALSE,col.names=NA)
  
  data$group<-rownames(data)
  data$group2<-data$group
  
  pdf(opts$outfile,width=3*9,height=6.5*n)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(n,3)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  for (i in 1:n) {
    gnum <- nlevels(as.factor(meta[,i]))
    #gcol <- rainbow(gnum)
	if(gnum < 22){
		gcol <- c(brewer.pal(3, "Set3")[1],brewer.pal(12, "Set3")[12],brewer.pal(12, "Set3")[3:11],brewer.pal(4, "Set3")[3],brewer.pal(3, "Pastel1")[1:2],brewer.pal(3,"Accent")[1],brewer.pal(8,"Set2")[6:8],brewer.pal(3,"Dark2"),length(gnum))
	}else{
		lc <- c(brewer.pal(3, "Set3")[1],brewer.pal(12, "Set3")[12],brewer.pal(12, "Set3")[3:11],brewer.pal(4, "Set3")[3],brewer.pal(3, "Pastel1")[1:2],brewer.pal(3,"Accent")[1],brewer.pal(8,"Set2")[6:8],brewer.pal(3,"Dark2"),length(21))
		gcol <- colorRampPalette(lc)(gnum)
	}
    #gshape <- sample(0:19,gnum)  
    data$group <- meta[,i] # grouping information
    
    if (opts$drawlabel) {
      umap12<-ggplot() +
        geom_point(data=data,aes(x=UMAP1,y=UMAP2,color=group),stat='identity',size=ps) + 
        geom_text(data=data,aes(x=UMAP1,y=UMAP2,label=group2),hjust=0,vjust=-1,alpha=0.8) +
        scale_colour_manual(values=gcol)+guides(col=guide_legend(colnames(meta)[i],ncol=ceiling(gnum/14))) +
        xlab("UMAP1") + ylab("UMAP2")  +
        theme(axis.text.x=element_text(size=14,colour="black"),
              axis.text.y=element_text(size=14,colour="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              panel.border=element_rect(fill=NA),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
			  legend.key = element_rect(fill="white"),
			  legend.text = element_text(size=14, color="black"),
			  legend.key.size = unit(0.8,"cm"),
			  legend.title = element_text(size = 16),
			  plot.margin = unit(rep(1.5,4),'lines'))
      
      
      umap13<-ggplot() +
        geom_point(data=data,aes(x=UMAP1,y=UMAP3,colour=group),stat='identity',size=ps) +
        geom_text(data=data,aes(x=UMAP1,y=UMAP3,label=group2),hjust=0,vjust=-1,alpha=0.8) +
        scale_colour_manual(values=gcol)+guides(col=guide_legend(colnames(meta)[i],ncol=ceiling(gnum/14))) +
        xlab("UMAP1") + ylab("UMAP3") + 
        theme(axis.text.x=element_text(size=14,colour="black"),
              axis.text.y=element_text(size=14,colour="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              panel.border=element_rect(fill=NA),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
			  legend.key = element_rect(fill="white"),
			  legend.text = element_text(size=14, color="black"),
			  legend.key.size = unit(0.8,"cm"),
			  legend.title = element_text(size = 16),
			  plot.margin = unit(rep(1.5,4),'lines'))
      
      umap23<-ggplot() +
        geom_point(data=data,aes(x=UMAP2,y=UMAP3,colour=group),stat='identity',size=ps) +
        geom_text(data=data,aes(x=UMAP2,y=UMAP3,label=group2),hjust=0,vjust=-1,alpha=0.8) +
        scale_colour_manual(values=gcol)+guides(col=guide_legend(colnames(meta)[i],ncol=ceiling(gnum/14))) +
        xlab("UMAP2") + ylab("UMAP3") + 
        theme(axis.text.x=element_text(size=14,colour="black"),
              axis.text.y=element_text(size=14,colour="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              panel.border=element_rect(fill=NA),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
			  legend.key = element_rect(fill="white"),
			  legend.text = element_text(size=14, color="black"),
			  legend.key.size = unit(0.8,"cm"),
			  legend.title = element_text(size = 16),
			  plot.margin = unit(rep(1.5,4),'lines'))
    } else {
      umap12<-ggplot() +
        geom_point(data=data,aes(x=UMAP1,y=UMAP2,colour=group),stat='identity',size=ps) +
        scale_colour_manual(values=gcol)+guides(col=guide_legend(colnames(meta)[i],ncol=ceiling(gnum/14))) +
        xlab("UMAP1") + ylab("UMAP2") + 
        theme(axis.text.x=element_text(size=14,colour="black"),
              axis.text.y=element_text(size=14,colour="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              panel.border=element_rect(fill=NA),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
			  legend.key = element_rect(fill="white"),
			  legend.text = element_text(size=14, color="black"),
			  legend.title = element_text(size = 16),
			  legend.key.size = unit(0.8,"cm"),
			  plot.margin = unit(rep(1.5,4),'lines'))
      
      
      umap13<-ggplot() +
        geom_point(data=data,aes(x=UMAP1,y=UMAP3,colour=group),stat='identity',size=ps) +
        scale_colour_manual(values=gcol)+guides(col=guide_legend(colnames(meta)[i],ncol=ceiling(gnum/14))) +
        xlab("UMAP1") + ylab("UMAP3") + 
        theme(axis.text.x=element_text(size=14,colour="black"),
              axis.text.y=element_text(size=14,colour="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              panel.border=element_rect(fill=NA),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
			  legend.key = element_rect(fill="white"),
			  legend.text = element_text(size=14, color="black"),
			  legend.key.size = unit(0.8,"cm"),
			  legend.title = element_text(size = 16),
			  plot.margin = unit(rep(1.5,4),'lines'))
      
      umap23<-ggplot() +
        geom_point(data=data,aes(x=UMAP2,y=UMAP3,colour=group),stat='identity',size=ps) +
        scale_colour_manual(values=gcol)+guides(col=guide_legend(colnames(meta)[i],ncol=ceiling(gnum/14))) +
        xlab("UMAP2") + ylab("UMAP3") +
        theme(axis.text.x=element_text(size=14,colour="black"),
              axis.text.y=element_text(size=14,colour="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              panel.border=element_rect(fill=NA),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
			  legend.key = element_rect(fill="white"),
			  legend.text = element_text(size=14, color="black"),
			  legend.key.size = unit(0.8,"cm"),
			  legend.title = element_text(size = 16),
			  plot.margin = unit(rep(1.5,4),'lines'))
    }
    
    
    print(umap12,vp=vplayout(i,1))
    print(umap13,vp=vplayout(i,2))
    print(umap23,vp=vplayout(i,3))
  }
  
  invisible(dev.off())
  
}

## calc
PM_umap(da=dst_orig, md=meta_orig)
