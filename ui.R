library(shiny)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

#######################raise max file upload limit#####################################

cnvFiles = list.files("Segments/", pattern="*cnv_genomicSegment_rmGermline")
cnvLocations <- as.list(rep("Segments/",length(cnvFiles)))
for (cnv in 1:length(cnvFiles)){
  cnvLocations[cnv] <- paste(cnvLocations[cnv],cnvFiles[cnv], sep="")
}
cancerTypes <- c(rep("x",length(cnvFiles)))
for (cancer in 1:length(cnvFiles)){
  cancerTypes[cancer] <- substr(cnvFiles[cancer], 1, nchar(cnvFiles[cancer])-40)
}

cancer_names <- as.data.frame((fread("cancer_list.csv", sep = ",")))$full_name


mutationFiles = list.files("Mutations/", pattern="*_mutation")
mutationLocations <- as.list(rep("Mutations/",length(mutationFiles)))
for (mutation in 1:length(mutationFiles)){
  mutationLocations[mutation] <- paste(mutationLocations[mutation],mutationFiles[mutation], sep="")
}

genelocations <- as.data.frame(fread(paste0("Annotations/hg38_singleSymbol.bed"),sep="\t",stringsAsFactors = FALSE,col.names = c("symbol","chrom","start","end")))
geneOptions <- genelocations$symbol

###User interface###

fluidPage(
  
  titlePanel("CAIRN: Copy Alterations Intuitive-Rendering Navigator"),
  
  sidebarPanel(
    selectInput('cancer', 'Specify Cancer Type', cancer_names, selected ='Ovarian Cancer (OV) [TCGA]'),
    fileInput('custom_CNA', 'Or, upload custom CNA segment file (tab or comma delimited)', multiple = FALSE),
    helpText('This may be used to enable the analysis of subsets of tumor types or other custom or unpublished datasets.',
             'You may click below for an example showing the correct data format.  Currently using hg38/GRCh38 for coordinates.'),
    downloadLink("downloadDataExample", "Download example custom input data"),
    hr(),
    selectInput('interest_type', 'Type of query', choices = c("genes", "ends", "overlaps")),
    conditionalPanel(
      condition = "input.interest_type == 'genes'",
      helpText('All queries must be present on the same chromosome.',
               'Genes will query all segments which include any part of the indicated gene(s).',
               'Please select your gene(s)'),
      textInput('genelist', 'Select Gene', value = "BRCA1"),
      textInput('genemultilist', 'Force segments to also include genes (enter as: GENE1, GENE2):', value = "No additional genes"),
      textInput('generemovelist', 'Remove any segments containing gene', value = "None")),
    conditionalPanel(condition = "input.interest_type == 'ends'",
                     helpText('Please input the genomic region where you want the output segments to end on.')),
    conditionalPanel(condition = "input.interest_type == 'overlaps'",
                     helpText('Please input the genomic region where you want the output segments to completely overlap.')),
    conditionalPanel(
      condition = "input.interest_type == 'ends' || input.interest_type == 'overlaps'",
      textInput('chromosome_num_input','Chromosome', value = "chr17"),
      conditionalPanel(condition = "input.centromeric == false && input.telocentric == false",
      textInput('chromosome_start_input','Start coordinate (hg38/GRCh38)', value = "4000000"),
      textInput('chromosome_end_input','End coordinate (hg38/GRCh38)', value = "5000000")),
      checkboxInput('remove_regions', 'Remove any segments which overlap a specific region', value = FALSE)),
    
    conditionalPanel(
      condition = "input.interest_type == 'ends'",
      checkboxInput('telocentric', 'Display segments which start/end near telomeres', value = FALSE),
      checkboxInput('centromeric', 'Display segments which start/end near centromeres', value = FALSE)), 

    conditionalPanel(
          condition = "input.remove_regions == true",
          checkboxInput('remove_centomeres', 'Remove any segments which cross centromere', value = FALSE),
          checkboxInput('remove_custom', 'Remove any segments which cross a custom region', value = FALSE)),
    conditionalPanel(    
      condition = "input.remove_custom == true",
    textInput('chromosome_start_input_remove','Start coordinate (hg38/GRCh38)', value = "6000000"),
    textInput('chromosome_end_input_remove','End coordinate (hg38/GRCh38)', value = "7000000")),
    
    tags$head(
      tags$style(HTML('#run{background-color:#9999FF; text-color:black; border-color:black; font-weight:bold}'))
    ),
    actionButton('run', 'Run CAIRN'),
          
    hr(),
    helpText('Other options'),
    selectInput('cutoff', 'Minimum CNA amplitude', choices = (1:10)/10, selected = 0.2),
    textInput('lengthCutoff', 'Minimum CNA length (bp)', value = "10000"),
    textInput('flanking_length', 'Extend query region endpoints (bp)', value = "1000"),

    checkboxInput('label_TSGs', 'Label COSMIC tumor suppressor genes', value = FALSE),
    checkboxInput('label_OGs', 'Label COSMIC oncogenes', value = FALSE),
    checkboxInput('mark_TSG_mutations', 'Mark mutant tumor suppressor genes on segments', value = FALSE),
    checkboxInput('mark_OG_mutations', 'Mark mutant oncogenes on segments', value = FALSE),
    conditionalPanel(    
      condition = "input.amps_only == false",
    checkboxInput('deletions_only', 'Plot only deletions', value = FALSE)),
    conditionalPanel(    
      condition = "input.deletions_only == false",
    checkboxInput('amps_only', 'Plot only amplifications', value = FALSE)),
    checkboxInput('constrain_plot', 'Specify exact plot output size', value = FALSE),
    conditionalPanel(    
      condition = "input.constrain_plot == true",
      textInput('x_pixels', 'Number of pixels (x)', value = 800),
      textInput('y_pixels', 'Number of pixels (y)', value = 800))

  ) ,
  
  mainPanel(
    helpText('This visualization tool explores copy number alterations discovered in published cancer datasets.'),
    helpText('It is intended to help the oncology community observe of the relative rates of amplification, deletion, and mutation of interesting genes and regions.'),
    helpText('For help finding the name of your gene, please visit: https://genome.ucsc.edu/'),
    helpText('Blue segments indicate deletions, red segments indicate amplifications.  Darker shades indicate stronger magnitude of alteration.'),
    textOutput('header'),
    hr(),
    textOutput('line1'),
    textOutput('line2'),
    textOutput('line3'),
    conditionalPanel(
      condition = "output.executed == 'done' ",
    downloadButton("downloadData", "Download Segments")),
    br(),
    uiOutput('ui_plot'),
    br(),
    br(),
    img(src = "CAIRNlegend.jpg", height = 200, width = 350)
    
  )
)


