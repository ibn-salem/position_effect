#########################################################################
# A script to reformat or join two imput tables with HPO terms linked 
# to patients and breakpoint coordinates linked to patients.
# It outputs files with exact breakpoint locations and windows around the
# breakpoints as well as control reginos, that are placed randomly in the genome.
#
# USAGE:
# 
# Rscript combine_breakpoints_and_phenotypes.R \
#   SUBJECT_TO_PHENOTYPES \
#   BREAKPOINTS_BED \
#   WINDOE_SIZE \
#   OUTPUT_FILE
#
# 2017 by Jonas Ibn-Salem <j.ibn-salem@uni-mainz.de>
#########################################################################

require(stringr) # for str_split_fixed()
require(readr)  # for tsv file writing

#-------------------------------------------------------------------------------
# Commandline arguments and parameters: 
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# input file with mapping of patient ID to phenotypes
# inFile <- "data/HPO_distal_cases.txt"
inFile <- args[1]

# input BED file with breakpoint coordinates
# inFileBreakPoints <- "non_coding_DGAP_positions.bed"
inFileBreakPoints <-  args[2]

# Window sizes in Mb
WIN_SIZE <- parse_integer(args[3])

# output file
outFileWindow <- args[4]

#-------------------------------------------------------------------------------
# parse input data
#-------------------------------------------------------------------------------

patientDF <- read.table(inFile, header = TRUE)

bpDF <- read.table(inFileBreakPoints)

#-------------------------------------------------------------------------------
# Combine breakpoint locations with phenotype data
#-------------------------------------------------------------------------------

names(bpDF) <- c("chr", "start", "end", "Breakpoint_ID")

# convert breakpoint ID to character vecotr
bpDF$Breakpoint_ID <- as.character(bpDF$Breakpoint_ID)

# extract patient IDs from Breakpoint IDs by using the suffix before underscore '_'
patientIDs <- str_split_fixed(bpDF$Breakpoint_ID, '_', 2)[,1]

# map patient IDs to HPO term and concatenate all terms for each patient with ';'inFile
bpDF$HPO <-  sapply(patientIDs, function(ID){
  paste(patientDF[patientDF$ID == ID,  "HPO"], collapse=";")
})

# reorder colums to match the input format of topodombar tool:
bpDF <- bpDF[,c("chr", "start", "end", "Breakpoint_ID", "HPO")]

# # write output file
# write.table(bpDF, file = outFileBP, col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)


#-------------------------------------------------------------------------------
# take predefined window around breakpoints
#-------------------------------------------------------------------------------

bpWinDF <- bpDF

# define center of input regions and extend by half of the window size to both directions
center <- round( (bpWinDF$start + bpWinDF$end) / 2 )
bpWinDF$start <- center - WIN_SIZE/2
bpWinDF$end <- center + WIN_SIZE/2

# write output file
# write.table(bpWinDF, file = paste0(outFileWindow, ".", WIN_SIZE, "MB_win.bed"),
#             col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write_tsv(bpWinDF, outFileWindow,
          col_names = FALSE)

