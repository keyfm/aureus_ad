# install.packages(c("SnowballC","wordcloud","RColorBrewer","gridExtra","ggplot2"),lib="/home/fkey/Rlib/") # added "/home/fkey/Rlib/" to ~/.Renviron
# library("tm")
library("SnowballC")
library("RColorBrewer")
library("wordcloud")
library("gridExtra")
library("ggplot2")

#args = commandArgs(trailingOnly=TRUE)
# read samples.csv that includes Subject identifier
samples <- read.csv('samples.csv',stringsAsFactors=FALSE,header=TRUE)


## Hard coded Variables (eventually later passed as argument)
cutoffPercRead = 80
tax <- c("Staphylococcus_aureus","Staphylococcus","Staphylococcaceae","Bacillales","Bacilli","Firmicutes","Bacteria")
names(tax) <- c("species","genus","family","order","class","phylum","domain")

## load species data
d <- read.table('2-kraken2/sumKrakTest_species.txt')
subset <- d[ d[,3]==tax['species'] ,  ]
subset[,4] <- as.character(subset[,4])

## build output folder for pdf's 2-spades
system("mkdir -p pdf 3-samplesFiltPerSubject")


####################################################################################
############## Filter Samples for species-of-interest and sort by patient
####################################################################################
## Get IDs of samples that fullfill cutoff 
## > output one file per Subject ID that includes seq IDs which fullfill quality criteria (ie: %reads on species >= cutoffPercRead )

patLs <- list()

patSplCount <- numeric(length(unique(samples[,'Subject']))+1) # define vector for counting samples per subject. +1 for unknown category (sample could not be assigned...sanity check)
names(patSplCount) <- sort(c(unique(samples[,'Subject']),'unknown')) # -"-
patPosSplCount <- numeric(length(patSplCount)) # define vector to count samples per subject that fullfill %reads on target
names(patPosSplCount) <- names(patSplCount)
for (i in 1:nrow(subset)){ # loop over each species-of-interest report
    patID <- as.character(samples[ which(samples[,'Sample']==subset[i,4]) , 'Subject']) # get patient ID for sample
    if( length(patID) > 0 ){ # test sample has assigned subject id or not and count either way
        patSplCount[[ patID ]] = patSplCount[[ patID ]] + 1
    } else { # sample has no subject assigned. should not happen!
            patSplCount[[ 'unknown' ]] = patSplCount[[ 'unknown' ]] + 1
    }
    if(subset[i,1] >= cutoffPercRead){ # if sample has sufficient on-target reads
        if( length(patID) > 0 ){            
            patLs[[ patID ]] <- c(patLs[[ patID ]] , subset[i,4])
            patPosSplCount[[ patID ]] = patPosSplCount[[ patID ]] + 1
        } else { # sufficient data but sample has no subject assigned. should not happen!
            patPosSplCount[[ 'unknown' ]] = patPosSplCount[[ 'unknown' ]] + 1
        }
    }
}


# write sample IDs that had sufficient %reads on target into per Subject file. NOTE: subjects with 0 successful samples will not be made!
for( n in names(patLs) ){
      write.table(x=patLs[[n]],file=paste("3-samplesFiltPerSubject/subject",n,"_samples.txt",sep=""),quote=F,row.names=F,col.names=F)
}

# write dummy file with all patients that fullfill criteria
write.table(x=names(patLs),file="3-samplesFiltPerSubject/SubjectsPassKraken.txt",quote=F,row.names=F,col.names=F)

# print(paste("IMPORTANT: From a total of",nrow(subset),"samples positive for",tax['species'],ctr,"harbour >=",(100*cutoffPercRead),"% species-specific reads and are kept for assembly."))


##########################################################################################
############### Plot Various Summaries to Assess Quality
##########################################################################################


###############
##### Barplot %samples that fullfilled %reads on target cutoff
###############
percPosSpl <- patPosSplCount/patSplCount # get freq of successfull samples per subject

if(length(percPosSpl)/2<4){plot_width <- 4} else {plot_width <- length(percPosSpl)/2} # scale width of barplot
pdf("pdf/KrakenAssess_percSamplePerSubject.pdf",width=plot_width,height=4)
ggplot(data.frame(percPosSpl),aes(seq_along(percPosSpl),percPosSpl))+
 geom_bar(stat="identity",fill="seagreen2") +
 geom_text(aes(label=patSplCount), vjust=0) +
 scale_x_continuous( breaks=seq(1,length(percPosSpl),1), labels=names(percPosSpl)) +
 scale_y_continuous(minor_breaks = seq(0 , 1, 0.2)) +
 labs(title="%Samples w/ suffient %reads-on-target per Subject (N top of bar)", x="Subject", y = "%Samples usable for Assmbly")+
 theme( panel.grid.minor = element_blank())
dev.off()

###############
##### Histograms of %reads on target species and homo_sapiens
###############
pdf("pdf/KrakenAssess_percReadHist_wordcloud.pdf",width=4,height=4)

# hist target species
dd <- d[d[,3]==tax['species'] ,1]
qplot(dd,
      geom="histogram",
      binwidth = 1,  
      main = paste("Histogram %Reads Sum ",tax['species']), 
      fill = "red",
      xlab = "%Reads on tax node",
      xlim = c(0,100)) + guides(fill=FALSE)

# hist homo (hard coded to screen for human contamination)
dd <- d[d[,3]=="Homo_sapiens" ,1]
qplot(dd,
      geom="histogram",
      binwidth = 1,  
      main = "Histogram %Reads Sum Homo sapiens", 
      xlab = "%Reads on tax node",
      xlim = c(-1,100)) 

########
##### Wordcloud of all species that also fullfill cutoff per sample
########
# wordcloud, note: highly abundant species excluded by code to maintain visibility in wordcloud
dd <- table(as.character(d[ d[,1]>cutoffPercRead ,3]))
wordcloud(words = names(dd), freq = dd, min.freq = 1,max.words=200, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()


###############
##### Scatter plot species vs. higher order tax levels
###############
# read data across all tax levels
ls <- list()
for (i in names(tax)){
    ls[[i]] <- read.table(paste('2-kraken2/sumKrakTest_',i,'.txt',sep=""),stringsAsFactors=F)
}

# s <- "Staphylococcus_aureus"
# g <- "Staphylococcus"
# f <- "Staphylococcaceae"
# o <- "Bacillales"
# c <- "Bacilli"
# p <- "Firmicutes"
# d <- "Bacteria"

# base dataset we want to compare the others to
sa <- ls[['species']][ ls[['species']][,3]==tax['species'] , ]
rownames(sa) <- sa[,4]


toiLs <- list()
for (t in names(tax)[-1]){ # -1 exclude species
    subLs <- ls[[ t ]][ ls[[ t ]][,3]==tax[ t ] , ] # get data only for tax ID of interest
    rownames(subLs) <- subLs[,4]
    subLsData <- vector()
    for (i in rownames(sa)){
        if( i %in% rownames(subLs) ){
            subLsData <- c(subLsData, subLs[i,1] ) # tax table and ID present, get %read value
        } else {
            subLsData <- c(subLsData, NA ) # tax table not ID present, get NA
        }
    }
    d <- data.frame(species=sa[,1],test=subLsData)
    toiLs[[t]] <- ggplot(d, aes(x = species, y = test)) + 
        geom_point()  + 
        ggtitle(paste("Plot ",t," vs. species (",tax[t]," vs. ",tax[1],")",sep="")) + 
         xlab(paste("%reads assigned",tax[1],sep=" ")) +
         ylab(paste("%reads assigned",tax[t],sep=" ")) +
	 theme_grey(base_size = 16)
}

pdf("pdf/KrakenAssess_percRead_Species_vs_varTaxNodes.pdf",width=7,height=7)
toiLs
dev.off()
