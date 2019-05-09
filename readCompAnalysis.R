#!/usr/bin/env Rscript

packages<-c("tidyverse","grid","gridExtra")
tmp<-lapply(packages,function(x) suppressPackageStartupMessages(require(x,character.only=T)))

Preprocess<-function(compFile,percentFilter=0.005)
	{
	##this function stores nucleotide composition per base data from a file in a data frame
	##subsets the data frame by a user given threshold to remove long tail, which is often extremely noisy
	##adds AT, GC columns to it
	# maxBases<-read_delim(compFile,"\t",skip=3,n_max=1)[,2]
	maxBases<-as.numeric(read_delim(compFile,"\t",skip=1,n_max=1)[,2])
	numFilter<-maxBases*percentFilter #number of base positions to subset until, avoids rendering noise at the end
	##**future enhancement: read the file in chunks and only do it for the data lines you need
	# df<-read_delim(compFile,"\t",skip=3)
	# colnames(df)<-c("POS","bases","A","C","G","T","N","avgQ","errQ","low","high")
	df<-read_delim(compFile,"\t")
	##**optimize this part by scanning Bases from the 1st data line and reading file into df till the limit specified
	df<-dplyr::filter(df,bases>=numFilter) #remove low confidence compositions by filtering the long tail
	df<-mutate(df,AT=A+T,GC=G+C)
	}

AreaPlot<-function(df,zoomLen=NULL)
	{
	##this function makes an area plot for given compositions per position
	##adds average base quality and sequence length distribution for full length plots
	##subsets data frame for zoomed plots
	if(!is.null(zoomLen)){df<-dplyr::filter(df,POS<=zoomLen)} #subsetting condition to follow for zoom plots
	ap<-ggplot(df,aes(x=POS,y=composition,fill=nt))+
		geom_area()+ #area plot
		theme_bw()+ #remove grey background
		geom_line(aes(y=avgQ),color="black",alpha=0.8) #avg Phred score
	# ap<-ap+theme(title=element_blank(),legend.position="none") #remove all title elements and legend
	if(is.null(zoomLen)) #condition to skip for zoom plots
		{
		ap<-ap+geom_line(aes(y=bases/max(bases)*100),color="red")+ #make normalized length freq as line
			scale_y_continuous(sec.axis=sec_axis(~.*max(df$bases)/100,name="Frequency")) #rev-norm freq for axis
		}
	ap<-ap+theme(title=element_blank(),legend.position="none") #remove all title elements and legend
	return(ap)
	}

###User defined parameters
zoom5Len<-300
zoom3Len<-100
fastqFile<-"ds2.fastq.gz" #fastq file to compute compositions on
rawQualFile<-"ds2-qual.txt" #file with normal base compositions; 5' aligned
revQualFile<-"ds2Compl-qual.txt" #file with reverse complement base compositions; 3' aligned
percentFilter<-0.005 #percent bases to filter until; parameter that avoids rendering noise at the end; chosen by the user
outFile<-"ds2-fastqCompositions.png" #file to export graphs in

###Process regular compositions
##Process full base compositons for center plots
rawComp<-Preprocess(rawQualFile,percentFilter) #sending to function for file reading
rawComp4nt<-gather(rawComp,key=nt,value=composition,c(A,T,G,C,N)) #gathers, reorders nt in long format to use with ggplot
rawComp4nt$nt<-factor(rawComp4nt$nt,levels=c("A","G","C","T","N")) #setting the order for nucleotide cols
areaRawComp4nt<-AreaPlot(rawComp4nt) #plotting full area graph with 4 nt
rawComp2nt<-gather(rawComp,key=nt,value=composition,c(AT,GC,N)) #reorders AT, GC in long format for ggplot
areaRawComp2nt<-AreaPlot(rawComp2nt) #plotting full area graph with combination compositions
##Process 5' end compositions for left zoomed-in plots
areaRawComp4ntZoom5<-AreaPlot(rawComp4nt,zoom5Len)
areaRawComp2ntZoom5<-AreaPlot(rawComp2nt,zoom5Len)

###Process reverse complement compositions
revComp<-Preprocess(revQualFile,percentFilter) #reading and filtering reverse complement compositions
revComp4nt<-gather(revComp,key=nt,value=composition,c(A,T,G,C,N)) #gathers, reorders nt in long format to use with ggplot
revComp4nt$nt<-factor(revComp4nt$nt,levels=c("A","G","C","T","N")) #setting the order for nucleotide cols
areaRevComp4ntZoom3<-AreaPlot(revComp4nt,zoom3Len)
areaRevComp4ntZoom3<-areaRevComp4ntZoom3+scale_x_reverse()
revComp2nt<-gather(revComp,key=nt,value=composition,c(AT,GC,N)) #reorders AT, GC in long format for ggplot
areaRevComp2ntZoom3<-AreaPlot(revComp2nt,zoom3Len)
areaRevComp2ntZoom3<-areaRevComp2ntZoom3+scale_x_reverse()

###Setting graphs in figure
areaArranged<-grid.arrange(areaRawComp4ntZoom5,areaRawComp4nt,areaRevComp4ntZoom3,areaRawComp2ntZoom5,areaRawComp2nt,areaRevComp2ntZoom3,nrow=2,ncol=3)
# centerPlots<-grid.arrange(areaRawComp4nt,areaRawComp2nt,ncol=1,top=textGrob(label="Aligned 5', length 1-XX nt"))
# leftPlots<-grid.arrange(areaRawComp4ntZoom5,areaRawComp2ntZoom5,ncol=1,top=textGrob(label="Aligned 5\' end"))
# rightPlots<-grid.arrange(areaRevComp4ntZoom3,areaRevComp2ntZoom3,ncol=1,top=textGrob(label="Aligned 3\' end"))

ggsave(outFile,plot=areaArranged)
