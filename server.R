library(shiny)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

options(shiny.maxRequestSize=300*1024^2)

###############NOTE: IF RUNNING LOCALLY, PLEASE RUN APP IN EXTERNAL MODE TO ENABLE DOWNLOADING OF DATA#############
#setwd("<path to folder where you downloaded CAIRN>")  ######set this too, then you're good to go!


#testing script, to enable editing of software line-by-line
#input <- list(interest_type = "ends", genelist = "BRCA1", genemultilist = "No additional genes",
#              chromosome_num_input = "chr17", chromosome_start_input = 3000000, chromosome_end_input = 4000000,
#              cutoff = 0.2, lengthCutoff = 10000, flanking_length = 1000,
#              generemovelist = "None", 
#              remove_centomeres = FALSE, label_TSGs = FALSE, label_OGs = FALSE,
#              mark_TSG_mutations = FALSE, mark_OG_mutations = FALSE,
#              cancer = "Ovarian Cancer (OV) [TCGA]",
#              telocentric = TRUE, centromeric = FALSE, remove_custom = FALSE)
#

#load constant data
genelocations <- as.data.frame(fread(paste0("Annotations/hg38_singleSymbol.bed"),sep="\t",stringsAsFactors = FALSE,col.names = c("symbol","chrom","start","end")))
species <- "human"
assembly <- "hg38"

#load constant functions

function(input, output) {
  values <- reactiveValues()
  
  example_download <- as.data.frame(fread(paste0("Segments/Example_segments"),col.names = c("sample","chrom","start","end","value")))
  output$downloadDataExample <- downloadHandler(
    filename = function() {
      paste("CNA_example.csv", sep="")
    },
    content = function(file) {
      write.csv(example_download, file,row.names=FALSE)
    },
    contentType = "text/csv"
  )
  
  CAIRNresult <- eventReactive(input$run,{
    
    #load NULL blank variables
    tidy_mut <- NULL
    colorval <- NULL
    mutations <- NULL
    cancerLabels <- as.data.frame((fread("cancer_list.csv", sep = ",")))
    
    if(species=="human"){
      chromEnds <- fread(paste0("Annotations/hg38_chromosomeEnds.bed"),stringsAsFactors = FALSE,col.names = c("chrom","start","end"))
      chromNames <- unique(chromEnds$chrom)
      centromeres <- as.data.frame(fread(paste0("Annotations/hg38_centromeres.bed"),stringsAsFactors = FALSE,col.names = c("chrom","start","end")))
    }
    if(species=="mouse"){
      chromEnds <- as.data.frame(fread(paste0("Annotations/mm10_chromosomeEnds.bed"),stringsAsFactors = FALSE,col.names = c("chrom","start","end")))
      chromNames <- unique(chromEnds$chrom)
      centromeres <- chromEnds}
    
    #load inputs
    if(is.null(input$custom_CNA)==FALSE){
      tumor_type <- input$custom_CNA$name
      cancerType <- input$custom_CNA$name
      if(input$mark_TSG_mutations ==TRUE | input$mark_OG_mutations ==TRUE){stop("Unable to graph mutation data with custom uploads")}
    }else{
      tumor_type <- cancerLabels$cancer[match(input$cancer, cancerLabels$full_name)]
      cancerType <- cancerLabels$cancer[match(input$cancer, cancerLabels$full_name)]}
    multiple_overlap <- FALSE
    
    if(input$interest_type == "genes"){
      genelist <- toupper(input$genelist)
      if(input$genemultilist != "No additional genes"){genelist <- toupper(unlist(c(genelist, strsplit(input$genemultilist, c(", ", " ",",")))))
      genelist <- unique(genelist[(genelist %in% genelocations$symbol)])
      multiple_overlap <- TRUE}
    }
    
    region_remove_list <- NULL
    if(input$generemovelist == "None"){generemovelist <- NULL} else {
      generemovelist <- input$generemovelist
    }
    
    if(input$label_TSGs == TRUE){label_TSGs <- TRUE} else {label_TSGs <- FALSE}
    if(input$label_OGs == TRUE){label_OGs <- TRUE} else {label_OGs <- FALSE}
    if(input$mark_TSG_mutations == TRUE){mark_TSG_mutations <- TRUE} else {mark_TSG_mutations <- FALSE}
    if(input$mark_OG_mutations == TRUE){mark_OG_mutations <- TRUE} else {mark_OG_mutations <- FALSE}
    
    inputSegmentsPath <- paste0("Segments/",tumor_type,"_TCGA_SNP6_cnv_genomicSegment_rmGermline")
    
    
    ############################################################################ 
    
    #load gene sites
    if(any(label_TSGs, label_OGs, mark_TSG_mutations, mark_OG_mutations)){
      genelocations <- as.data.frame(fread(paste0("Annotations/hg38_singleSymbol.bed"),sep="\t",stringsAsFactors = FALSE,col.names = c("symbol","chrom","start","end")))
      
      #load tumor suppressor genes (TSGs) and oncogenes (OGs)
      tumor_genes <- as.data.frame(fread("Annotations/COSMIC_TSG_OG.csv",stringsAsFactors = FALSE))
    }
    
    #load mutation files if any true
    if(any(mark_TSG_mutations, mark_OG_mutations)){
      mutations <- as.data.frame(fread(paste0("Mutations/",tumor_type,"_mutation"),stringsAsFactors = FALSE))
      colnames(mutations) <- substr(colnames(mutations),1,12)
      mutations <- subset(mutations, sample %in% tumor_genes$symbol)
      genes_to_show_mutation <- unique(c(if(mark_TSG_mutations==TRUE){subset(tumor_genes, tsg_og == "TSG")$symbol},if(mark_OG_mutations==TRUE){subset(tumor_genes, tsg_og == "OG")$symbol}))
      samples_with_mut_data <- colnames(mutations)[2:ncol(mutations)]
    }
    
    #user data: load interest sites files
    if(input$interest_type == "ends"){
      
      interestSites <- as.data.frame(matrix(0,nrow=1,ncol=3))
      colnames(interestSites) <- c("chrom", "start", "end")
      if(class(input$chromosome_num_input)=="integer"){
        interestSites$chrom <- paste0("chr",input$chromosome_num_input)}else{
          interestSites$chrom <- tolower(input$chromosome_num_input)
        }
      if(tolower(input$chromosome_num_input) == "chrx"){interestSites$chrom <- "chrX"}
      if(tolower(input$chromosome_num_input) == "chry"){interestSites$chrom <- "chrY"}
      selected_chrom <- interestSites$chrom
      if(input$telocentric == TRUE){
        interestSites <- subset(as.data.frame(chromEnds),chrom==selected_chrom)
        interestSites[2,] <- interestSites[1,]
        interestSites$start[2] <- interestSites$end[1] - 2000000 - as.integer(input$flanking_length)
        interestSites$end[1] <- interestSites$start[1] + 2000000 + as.integer(input$flanking_length)
      }
      if(input$centromeric == TRUE){
        if(input$telocentric == TRUE){teloSites <- interestSites}
        interestSites <- subset(as.data.frame(centromeres),chrom==selected_chrom)
        interestSites$end <- interestSites$end +2000000 + as.integer(input$flanking_length)
        interestSites$start <- interestSites$start - 2000000 - as.integer(input$flanking_length)
        if(input$telocentric == TRUE){interestSites <- rbind(teloSites,interestSites)}
      }
      
      
      if(sum(input$telocentric, input$centromeric)==0){
        interestSites$start <- as.integer(input$chromosome_start_input) - as.integer(input$flanking_length)
        interestSites$end <- as.integer(input$chromosome_end_input) + as.integer(input$flanking_length)}
    }
    
    if(input$interest_type == "overlaps"){
      
      interestSites <- as.data.frame(matrix(0,nrow=1,ncol=3))
      colnames(interestSites) <- c("chrom", "start", "end")
      if(class(input$chromosome_num_input)=="integer"){
        interestSites$chrom <- paste0("chr",input$chromosome_num_input)}else{
          interestSites$chrom <- tolower(input$chromosome_num_input)
        }
      if(tolower(input$chromosome_num_input) == "chrx"){interestSites$chrom <- "chrX"}
      if(tolower(input$chromosome_num_input) == "chry"){interestSites$chrom <- "chrY"}
      interestSites$start <- as.integer(input$chromosome_start_input) - as.integer(input$flanking_length)
      interestSites$end <- as.integer(input$chromosome_end_input) + as.integer(input$flanking_length)
    }
    
    if(input$remove_custom==TRUE){
      region_remove_list <- as.data.frame(matrix(0,nrow=1,ncol=3))
      colnames(region_remove_list) <- c("chrom", "start", "end")
      region_remove_list$chrom <- interestSites$chrom
      region_remove_list$start <- as.integer(input$chromosome_start_input_remove)
      region_remove_list$end <- as.integer(input$chromosome_end_input_remove)
    }
    
    if(input$interest_type == "genes"){
      if(length(genelist)==1){multiple_overlap <- FALSE}
      interestSites <- genelocations
      if(length(generemovelist)>0){removalSites <- subset(interestSites, symbol %in% generemovelist)
      removalSites <- removalSites[,c("chrom","start","end")]}
      interestSites <- subset(interestSites, symbol %in% genelist)
      if(multiple_overlap == TRUE){
        if(length(unique(interestSites$chrom)) == 1){
          interestSites[1,] <- c("all_genes",interestSites$chrom[1], min(interestSites$start, interestSites$end), max(interestSites$start, interestSites$end), "-",1)
          interestSites <- interestSites[1,]
        } else {stop("genes are not all on the same chromosome")
          interestSites <- NULL}}
      graphLabels <- interestSites$symbol
      genelist <- graphLabels
      interestSites <- interestSites[,c("chrom","start","end")]
      interestSites$start <- as.integer(interestSites$start)
      interestSites$end <- as.integer(interestSites$end)
    }
    
    #user data: load segments file
    if(is.null(input$custom_CNA)==FALSE){segments <- as.data.frame(fread(input$custom_CNA$datapath,stringsAsFactors = FALSE))}else{
      segments <- as.data.frame(fread(inputSegmentsPath,stringsAsFactors = FALSE))}
    if(ncol(segments)==3){
      colnames(segments) <- c("chrom","start","end")
      segments$sample <- "sample1"
      segments$value <- 1
      segments <- segments[,c("sample","chrom","start","end", "value")]
    } else {if(ncol(segments)>4){
      colnames(segments) <- c("sample","chrom","start","end","value")
      segments <- subset(segments, abs(value)>input$cutoff)
      values_exist <- TRUE
    } else {colnames(segments) <- c("sample","chrom","start","end")}}
    if(class(segments$chrom)=="integer"){segments$chrom <- paste0("chr",segments$chrom)}
    if(species=="human"){segments$chrom <- gsub("23","X",segments$chrom)
    segments$chrom <- gsub("24","Y",segments$chrom)
    if(input$remove_centomeres == TRUE){region_remove_list <- rbind(centromeres,region_remove_list)} 
    if(sum(grepl("chr",segments$chrom))==0){segments$chrom <- paste0("chr",segments$chrom)}}
    if(species=="mouse"){segments$chrom <- gsub("20","X",segments$chrom)
    segments$chrom <- gsub("21","Y",segments$chrom)
    #chromEnds <- as.data.frame(fread(paste0("Annotations/mm10_chromosomeEnds.bed"),stringsAsFactors = FALSE,col.names = c("chrom","start","end")))
    #chromNames <- unique(chromEnds$chrom)
    if(sum(grepl("chr",segments$chrom))==0){segments$chrom <- paste0("chr",segments$chrom)}}
    
    if(any(mark_TSG_mutations, mark_OG_mutations)){segments <- subset(segments, sample %in% samples_with_mut_data)}
    
    segments$length <- segments$end - segments$start
    segments <- subset(segments, length >= as.numeric(input$lengthCutoff))
    segments <- segments[,1:(ncol(segments)-1)]
    all_samples <- unique(segments$sample)
    
    #obtain and graph segments of interest
    site=1
    interest_site <- NULL
    
    CAIRNsegments <- lapply(1:nrow(interestSites), function(site){
      
      interest_site <- interestSites[site,]
      removalSitesChr <- NULL
      if(input$interest_type == "genes" & is.null(generemovelist) == FALSE){removalSitesChr <- subset(removalSites,chrom == interest_site$chrom)}
      if(input$interest_type == "genes" & is.null(generemovelist) == TRUE & is.null(region_remove_list)==FALSE){removalSitesChr <- subset(region_remove_list,chrom == interest_site$chrom) }
      if((input$interest_type == "ends" | input$interest_type == "overlaps") & is.null(region_remove_list) == FALSE){removalSitesChr <- subset(region_remove_list,chrom == interest_site$chrom) }
      
      if(input$interest_type == "ends"){
        #subset segments into those ending at interest site
        segments_of_interest <- subset(segments, chrom == interest_site$chrom)
        segments_of_interest$start_overlap <- as.integer((segments_of_interest$start < interest_site$end) & (segments_of_interest$start > interest_site$start))
        segments_of_interest$end_overlap <- as.integer((segments_of_interest$end < interest_site$end) & (segments_of_interest$end > interest_site$start))
        if(is.null(removalSitesChr)){}else{
          if(length(region_remove_list)>0 & nrow(removalSitesChr) > 0){
            segments_of_interest$removeoverlap <- 0
            for(removal_site in 1:nrow(removalSitesChr)){
              remove_site <- removalSitesChr[removal_site,]
              segments_of_interest$removeoverlap <- segments_of_interest$removeoverlap +
                as.integer((remove_site$start > segments_of_interest$start) & (remove_site$start < segments_of_interest$end)
                           | (remove_site$end < segments_of_interest$end) & (remove_site$end > segments_of_interest$start))
            }
            removedSegmentNum <- sum(segments_of_interest$removeoverlap >0)
            removedSegmentPercent <- removedSegmentNum / nrow(segments_of_interest)
            segments_of_interest <- subset(segments_of_interest, removeoverlap == 0)[,1:(ncol(segments_of_interest)-1)]
          }}
        segments_of_interest$length <- segments_of_interest$end - segments_of_interest$start
        segments_of_interest1 <- subset(segments_of_interest, end_overlap>0)
        segments_of_interest1 <- segments_of_interest1[order(segments_of_interest1$length, decreasing=TRUE),]
        segments_of_interest2 <- subset(segments_of_interest, start_overlap>0)
        segments_of_interest2 <- segments_of_interest2[order(segments_of_interest2$length, decreasing=TRUE),]
        segments_of_interest <- rbind(segments_of_interest1,segments_of_interest2)
        segments_of_interest <- segments_of_interest[order(segments_of_interest$length, decreasing=TRUE),]
        if(nrow(segments_of_interest)==0){}else{
          segments_of_interest$yval <- sequence(nrow(segments_of_interest))
          segments_of_interest <- segments_of_interest %>% rowwise() %>% mutate (colorval = if(value<0){"blue"}else{"red"})
        }}
      
      if(input$interest_type == "overlaps"){
        segments_of_interest <- subset(segments, chrom == interest_site$chrom)
        segments_of_interest$overlap <- as.integer((interest_site$start > segments_of_interest$start) & (interest_site$start < segments_of_interest$end)
                                                   & (interest_site$end < segments_of_interest$end) )
        segments_of_interest <- subset(segments_of_interest, overlap >0)
        if(nrow(segments_of_interest)==0){}else{
          if(is.null(removalSitesChr)){}else{
            if(length(region_remove_list)>0 & nrow(removalSitesChr) > 0){
              segments_of_interest$removeoverlap <- 0
              for(removal_site in 1:nrow(removalSitesChr)){
                remove_site <- removalSitesChr[removal_site,]
                segments_of_interest$removeoverlap <- segments_of_interest$removeoverlap +
                  as.integer((remove_site$start > segments_of_interest$start) & (remove_site$start < segments_of_interest$end)
                             | (remove_site$end < segments_of_interest$end) & (remove_site$end > segments_of_interest$start))
              }
              removedSegmentNum <- sum(segments_of_interest$removeoverlap >0)
              removedSegmentPercent <- removedSegmentNum / nrow(segments_of_interest)
              segments_of_interest <- subset(segments_of_interest, removeoverlap == 0)[,1:(ncol(segments_of_interest)-1)]
            }}}
        if(nrow(segments_of_interest)==0){}else{
          segments_of_interest$length <- segments_of_interest$end - segments_of_interest$start
          segments_of_interest <- segments_of_interest[order(segments_of_interest$length, decreasing=TRUE),]
          segments_of_interest$yval <- sequence(nrow(segments_of_interest))
          segments_of_interest <- segments_of_interest %>% rowwise() %>% mutate (colorval = if(value<0){"blue"}else{"red"})
        }}
      
      if(input$interest_type == "genes"){
        segments_of_interest <- subset(segments, chrom == interest_site$chrom)
        segments_of_interest$overlap <- as.integer((interest_site$start > segments_of_interest$start) & (interest_site$start < segments_of_interest$end)
                                                   & (interest_site$end < segments_of_interest$end))
        segments_of_interest <- subset(segments_of_interest, overlap >0)
        if(nrow(segments_of_interest)==0){}else{
          segments_of_interest$length <- segments_of_interest$end - segments_of_interest$start
          segments_of_interest <- segments_of_interest[order(segments_of_interest$length, decreasing=TRUE),]
          if(is.null(generemovelist) & is.null(region_remove_list)){}else{
            if(is.null(generemovelist) & is.null(region_remove_list)==FALSE){
              if(length(region_remove_list)>0 & nrow(removalSitesChr) > 0){
                segments_of_interest$removeoverlap <- 0
                for(removal_site in 1:nrow(removalSitesChr)){
                  remove_site <- removalSitesChr[removal_site,]
                  segments_of_interest$removeoverlap <- segments_of_interest$removeoverlap +
                    as.integer((remove_site$start > segments_of_interest$start) & (remove_site$start < segments_of_interest$end)
                               | (remove_site$end < segments_of_interest$end) & (remove_site$end > segments_of_interest$start))
                }
                removedSegmentNum <- sum(segments_of_interest$removeoverlap >0)
                removedSegmentPercent <- removedSegmentNum / nrow(segments_of_interest)
                segments_of_interest <- subset(segments_of_interest, removeoverlap == 0)[,1:(ncol(segments_of_interest)-1)]
              }} else{
                if(length(generemovelist)>0 & nrow(removalSitesChr) > 0){
                  segments_of_interest$removeoverlap <- 0
                  for(removal_site in 1:nrow(removalSitesChr)){
                    remove_site <- removalSitesChr[removal_site,]
                    segments_of_interest$removeoverlap <- segments_of_interest$removeoverlap +
                      as.integer((remove_site$start > segments_of_interest$start) & (remove_site$start < segments_of_interest$end)
                                 | (remove_site$end < segments_of_interest$end) & (remove_site$end > segments_of_interest$start))
                  }
                  removedSegmentNum <- sum(segments_of_interest$removeoverlap >0)
                  removedSegmentPercent <- removedSegmentNum / nrow(segments_of_interest)
                  segments_of_interest <- subset(segments_of_interest, removeoverlap == 0)[,1:(ncol(segments_of_interest)-1)]
                  if(nrow(segments_of_interest)>1){}else{stop("no segments were found")}
                }}}
          segments_of_interest$yval <- sequence(nrow(segments_of_interest))
          segments_of_interest <- as.data.frame(segments_of_interest %>% rowwise() %>% mutate (colorval = if(value<0){"blue"}else{"red"}))
        }}
      
      if(nrow(segments_of_interest)>1){
        cut <-input$cutoff
        segments_of_interest <- segments_of_interest[((abs(segments_of_interest$value)) > 0.1),]
        segments_of_interest <- segments_of_interest[!duplicated(as.data.frame(segments_of_interest)[,1:5]),]
        segments_of_interest$yval <- sequence(nrow(segments_of_interest))#remove any duplicated segments
        segments_of_interest_export <- segments_of_interest
        segments_of_interest$value <-  (segments_of_interest$value > -1)*segments_of_interest$value + (segments_of_interest$value <= -1)*(-1) #set minimum CNA to -1 for graphing color consistency
        segments_of_interest$value <-  (segments_of_interest$value < 2)*segments_of_interest$value + (segments_of_interest$value >= 2)*(2)   #set maximum CNA to 2 for graphing color consistency    
      }
      return(list(segs = segments_of_interest, rem_segs = removalSitesChr))
    })
    
    segments_of_interest <- CAIRNsegments[[1]]$segs
    if(length(CAIRNsegments)>1){
      for(site in 2:nrow(interestSites)){
        segments_of_interest <- rbind(segments_of_interest, CAIRNsegments[[site]]$segs)
      }}
    
    removalSitesChr <- CAIRNsegments[[1]]$rem_segs
    if(length(CAIRNsegments)>1){
      for(site in 2:nrow(interestSites)){
        segments_of_interest <- rbind(segments_of_interest, CAIRNsegments[[site]]$rem_segs)
      }}
    
    
    if(nrow(segments_of_interest)>1){
      if(input$deletions_only==TRUE){segments_of_interest <- subset(as.data.frame(segments_of_interest), value<0)}
      if(input$amps_only==TRUE){segments_of_interest <- subset(as.data.frame(segments_of_interest), value>0)}
      
      sample_summary <- as.data.frame(matrix(0,nrow = length(unique(segments_of_interest$sample)), ncol=2))
      colnames(sample_summary) <- c("sample","length_segs")
      sample_summary$sample <- unique(segments_of_interest$sample)
      sample_summary$length_segs <- sapply(1:nrow(sample_summary), function(sample_num){
        return(sum(subset(segments_of_interest, sample == sample_summary$sample[sample_num])$length))
      })
      sample_summary <- sample_summary[order(sample_summary$length_segs, decreasing = TRUE),]
      sample_graph_order <- sample_summary$sample
      
      segments_of_interest$yval <- match(segments_of_interest$sample, sample_graph_order)
    }else{stop("no segments were found")}
    
    #generate labels for cancer genes
    if(any(label_TSGs, label_OGs)){
      tumor_genes$chrom <- genelocations$chrom[match(tumor_genes$symbol, genelocations$symbol)]
      chr_tumor_genes <- subset(tumor_genes, chrom == interestSites$chrom[1])
      chr_tumor_genes$start <- genelocations$start[match(chr_tumor_genes$symbol, genelocations$symbol)]
      chr_tumor_genes$end <- genelocations$end[match(chr_tumor_genes$symbol, genelocations$symbol)]
    }
    
    #generate mutation plotting data
    if(any(mark_TSG_mutations, mark_OG_mutations)){
      possible_mutations <- subset(genelocations, chrom %in% interestSites$chrom[1])
      possible_mutations <- subset(possible_mutations, symbol %in% genes_to_show_mutation)
      mutations <- subset(mutations, sample %in% possible_mutations$symbol)
      
      if(nrow(mutations)>0){
        tidy_mut <- gather(data= mutations, key="gene", value="mutant",  match(samples_with_mut_data, colnames(mutations)))
        colnames(tidy_mut) <- c("gene", "sample", "mutant")
        tidy_mut <- subset(as.data.frame(tidy_mut), mutant != 0)
        if(nrow(tidy_mut) > 0){
        tidy_mut$chrom <- interestSites$chrom[1]
        tidy_mut$start <- genelocations$start[match(tidy_mut$gene, genelocations$symbol)]
        tidy_mut$end <- genelocations$end[match(tidy_mut$gene, genelocations$symbol)]
        tidy_mut$yval <- segments_of_interest$yval[match(tidy_mut$sample, segments_of_interest$sample)]
        tidy_mut <- subset(tidy_mut, yval > 0)
        tidy_mut$overlap <- 0
        temp_seg <- segments_of_interest[1,]
        for (mutation in 1:nrow(tidy_mut)){ #check each potential mutation against segment changes, keep overlaps
          temp_seg <- subset(segments_of_interest, sample == tidy_mut$sample[mutation])
          if(temp_seg$start <= tidy_mut$start[mutation] & temp_seg$end >= tidy_mut$end[mutation]){tidy_mut$overlap[mutation] <- 1}
        }
        if(sum(tidy_mut$overlap)==0){stop("No tumor suppressor or oncogene mutations found in output segments")}
        tidy_mut <- subset(tidy_mut, overlap ==1 )
        
        tidy_mut$colorval <- "black"
        for (mutation in 1:nrow(tidy_mut)){
          if(tidy_mut$gene[mutation] %in% subset(tumor_genes, tsg_og == "OG")$symbol){tidy_mut$colorval[mutation] <- "green"}
        }
        } else {tidy_mut <- NULL}
      }}
    
    #generate chromosome graphic data
    interest_chromosome <- as.data.frame(subset(chromEnds, chrom== interestSites$chrom[1]))
    if(species == "human"){
      interest_centromere <- subset(centromeres, chrom== interestSites$chrom[1])
    }
    
    if(nrow(segments_of_interest)==0){}else{
      
      segPlot <- ggplot(segments_of_interest) +
        theme(plot.background = element_blank()
              ,panel.grid.major = element_blank()
              ,panel.grid.minor = element_blank()
              ,panel.border = element_blank()
              ,panel.background = element_blank()
              ,legend.position="none"
              ,axis.line= element_blank()
              ,axis.text= element_blank()
              ,axis.ticks= element_blank()
        )
      
      if(nrow(subset(segments_of_interest,colorval=="blue"))>0){
        segPlot <- segPlot + geom_segment(data = subset(segments_of_interest,colorval=="blue"),  aes(x=start, y=yval,xend = end, yend=yval, alpha=-value), size=.95, color="blue")
      }
      
      if(nrow(subset(segments_of_interest,colorval=="red"))>0){
        segPlot <- segPlot +   geom_segment(data = subset(segments_of_interest,colorval=="red"),  aes(x=start, y=yval,xend = end, yend=yval, alpha=value), size=.95, color="red")
      }
      
      segPlot <- segPlot +
        geom_segment(data = subset(chromEnds, chrom== interestSites$chrom[1]), aes(x=start, y=-1,xend = end, yend=-1), size=3 ,color="#666666")+
        geom_point(x=interest_chromosome$start, y=-1, size= 2.4,color="#999999")+
        geom_point(x=interest_chromosome$end, y=-1, size= 2.4,color="#999999")+
        labs(x="", y="Samples", title=if(input$interest_type=="genes"){paste0(genelist[site],"(",interestSites$chrom[1],":", min(interestSites$start),"-", max(interestSites$end),")")}else{paste0(interestSites$chrom[1],":", min(interestSites$start),"-", max(interestSites$end))})
      
      if(as.numeric(interestSites$end[1]) - as.numeric(interestSites$start[1]) > 1000000){
        segPlot <- segPlot +
          geom_segment(data = interestSites, aes(x=start, y=-1,xend =end, yend=-1), size=3 ,color="#009975")
      } else {
        segPlot <- segPlot + geom_point(x= interestSites$start, y =-1, size=2, color="#009975")
      }  
      
      if(is.null(tidy_mut)==FALSE){
        if(nrow(subset(tidy_mut, colorval == "green"))>0){ 
          segPlot <- segPlot + geom_point(data = subset(tidy_mut, colorval == "green"), aes(x= start, y = yval), color= "#4CBB17", size = 1.5, shape=17)
        }
        if(nrow(subset(tidy_mut, colorval == "black"))>0){ 
          segPlot <- segPlot + geom_point(data = subset(tidy_mut, colorval == "black"), aes(x= start, y = yval), color= "#000000", size = 1.5, shape=17)
        }}
      
      if(any(label_TSGs)){
        segPlot <- segPlot + geom_text(data = subset(chr_tumor_genes, tsg_og == "TSG"), aes(x= -1, y = -(10+max(segments_of_interest$yval)/30), label = symbol, angle=90), color= "#FFFFFF", size = 4)
        segPlot <- segPlot + geom_text(data = subset(chr_tumor_genes, tsg_og == "TSG"), aes(x= start, y = -(3+max(segments_of_interest$yval)/30), label = symbol, angle=90), color= "#000000", size = 4)
      }
      if(any(label_OGs)){
        segPlot <- segPlot + geom_text(data = subset(chr_tumor_genes, tsg_og == "OG"), aes(x= start, y = -(3+max(segments_of_interest$yval)/30), label = symbol, angle=90), color= "#4CBB17", size = 4)
      }
      
      
      if(species=="human"){
        segPlot <- segPlot + geom_segment(data = subset(centromeres, chrom== interestSites$chrom[1]), aes(x=start, y=-1,xend = end, yend=-1), size=3 ,color="#000000")
      }
      
      if(is.null(removalSitesChr)==FALSE){
        if(nrow(removalSitesChr)>0){
          for(removal_site in 1:nrow(removalSitesChr)){
            remove_site <- removalSitesChr[removal_site,]
            if(remove_site$end - remove_site$start > 1000000){
              segPlot <- segPlot +
                geom_segment(data = remove_site, aes(x=start, y=-1,xend =end, yend=-1), size=3 ,color="#FF0080")
            } else {
              segPlot <- segPlot + geom_point(x= remove_site$start, y =-1, size=2, color="#FF0080")
            }  
          }}}
    }
    
    
    ############################################################################
    if(is.null(mutations)==FALSE){
      sample_population <- unique(segments$sample, mutations$sample)}else{
        sample_population <- unique(segments$sample)
      }
    numDeletions <- sum(segments_of_interest$value < 0)
    numAmps <- sum(segments_of_interest$value > 0)
    percentSegDeletions <- round(numDeletions/ nrow(segments_of_interest) * 100, digits=1)
    percentSegAmps <- round(numAmps/ nrow(segments_of_interest) * 100, digits=1)
    percentSampleDeletions <- round(length(unique(subset(segments_of_interest, value<0)$sample)) / length(sample_population) * 100, digits=1)
    percentSampleAmps <- round(length(unique(subset(segments_of_interest, value>0)$sample)) / length(sample_population) * 100, digits=1)
    
    output$downloadData <- downloadHandler(
      filename = "CAIRNdownload.csv",
      content = function(file) {
        write.csv(segments_of_interest[,1:5], file, row.names=FALSE)
      },
      contentType = "text/csv"
    )

    
    
    return(list(cancer= cancerType, segPlot_vis = segPlot, numSeg = nrow(segments_of_interest),
                numTumors = length(sample_population), numDeletions = numDeletions, numAmps = numAmps,
                percentSegDeletions = percentSegDeletions, percentSegAmps = percentSegAmps, done = "done",
                percentSampleDeletions = percentSampleDeletions, percentSampleAmps = percentSampleAmps))
  })
  
  plot_height <- function() {
    # calculate values$facetCount
    values$segmentval <- CAIRNresult()$numSeg*1.5 + 400
    return(values$segmentval)
  }
  
  output$header <- renderText({
    paste(if(is.null(input$custom_CNA)){input$cancer}else{"Custom Data"}, "-", "Copy Number Alterations - Segment Analysis")
  })
  output$segPlot_vis <- renderPlot(CAIRNresult()$segPlot_vis)
  output$ui_plot <- renderUI({plotOutput("segPlot_vis", height = if(input$constrain_plot == TRUE){input$y_pixels}else{plot_height()}, width = if(input$constrain_plot == TRUE){input$x_pixels}else{"800px"})})
  output$line1 <- renderText(paste0(CAIRNresult()$numSeg, " CNAs were found in ", CAIRNresult()$numTumors, " samples with data are altered in query region."))
  output$line2 <- renderText(paste0(CAIRNresult()$numDeletions, " CNAs were deletions (", CAIRNresult()$percentSegDeletions, "% of CNAs, ",CAIRNresult()$percentSampleDeletions,"% of samples)"))
  output$line3 <- renderText(paste0(CAIRNresult()$numAmps, " CNAs were amplifications (", CAIRNresult()$percentSegAmps, "% of CNAs, ",CAIRNresult()$percentSampleAmps,"% of samples)"))
  
  output$executed <- renderText(CAIRNresult()$done)
  outputOptions(output, "executed", suspendWhenHidden=FALSE)
}
