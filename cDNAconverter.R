# Code to convert MANE genomic coordinates to complementary DNA coordinates.

# Necessary to calculate exon junction complex positions since MANE coordinates
# are not consecutive (i.e., implicitly includes introns). Also
# calculates instances of intron positions between stop codon and
# first 3' UTR exon (classified as "stop_intron").

# A full table of stop introns generated by this code can be found in
# the Data folder of this repository (3Prime_UTRs/Data/threeUTR_stopIntrons.tsv).

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(readr))
suppressMessages(library(dplyr))

# -------------------------------------------------------------------
# Download and read in Ensembl v110 .gtf file
# https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/
# Downloaded on August 1, 2023
ensemblGTF <- rtracklayer::import("C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/3C/3Prime_UTRs/Data/Input/Homo_sapiens.GRCh38.110.chr.gtf.gz")
ensemblDT  <- as.data.table(ensemblGTF)

# Pull all stop codons and 3' UTR exons
cdna_converter <- ensemblDT[type %in% c("stop_codon", "three_prime_utr") & tag == "MANE_Select", 
                            .(gene_name, transcript_id, type, seqnames, strand, start, end, exon_number, width)]

# Create "exon_number" values for 3' UTRs using stop_codon value
# (Exon numbers are not included in Ensembl .gtf file for UTRs)
fill_exonNumber <- function(dt) {
  stop_codon_exon <- as.numeric(dt[type == "stop_codon", exon_number])
  dt[type == "three_prime_utr", exon_number := seq.int(stop_codon_exon + 1, length.out = .N)]
  return(dt)
}

# Apply the function by grouping with gene_name and transcript_id
cdna_converter <- cdna_converter[, fill_exonNumber(copy(.SD)), by = .(gene_name, transcript_id)]

# Set the order by exon number (so we don't have to worry about strand)
cdna_converter[, exon_number := as.numeric(exon_number)]
setorder(cdna_converter, transcript_id, exon_number)

# Add the cumulative width of exons up to that exon number
cdna_converter <- cdna_converter[, .(cumwidth = cumsum(width), start = start, end = end,
                                     width = width, exon_number = exon_number),
                                by = .(gene_name, transcript_id, type, seqnames, strand)]

# Update the exonStart_cDNA and exonEnd_cDNA positions
cdna_converter[, exonStart_cDNA := cumwidth - width + 1]
cdna_converter[, exonEnd_cDNA := cumwidth]

# Update the exonStart_cDNA position for "three_prime_utr" rows
# Need exonStart_cDNA positions to begin from exonEnd_cDNA of stop_codon
for (i in 2:nrow(cdna_converter)) {
  message(sprintf("%i of %i", i, nrow(cdna_converter)))
  if(cdna_converter$type[i] == "three_prime_utr"){
    prev_exon_end <- cdna_converter$exonEnd_cDNA[i - 1]
    cdna_converter$exonStart_cDNA[i] <- prev_exon_end + 1
    cdna_converter$exonEnd_cDNA[i]   <- prev_exon_end + cdna_converter$width[i]
  }
}

# Output directory
outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/3C/3Prime_UTRs/Data"
outDir <- sprintf("%s/Data/cDNAconverter.R", dirname(outpath))
if(!file.exists(outDir)) dir.create(outDir)

write.table(cdna_converter, file = sprintf("%s/cdna_converter_ens110.tsv", outDir),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -------------------------------------------------------------------
# Import MANE 3' UTR exon and stop codon data
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/3C/3Prime_UTRs/Data/3primeUTR_NMD.R"
inDir  <- sprintf("%s/3primeUTR_NMD.R", dirname(inpath))

threeUTR_intronStop <- read_delim(sprintf("%s/threeUTR_intronStop.tsv", inDir), delim = "\t")
threeUTR_intronStop <- as.data.frame(threeUTR_intronStop)

# Merge input with cdna_converter.tsv file
# This takes care of stop codons/transcripts w/o 3' UTRs
selected_columns <- c("transcript_id", "start", "end", "phase")

convertedCDNA <- left_join(threeUTR_intronStop[selected_columns], cdna_converter, 
                           by = c("transcript_id", "start", "end"))

# Specify which coordinates are genomic coordinates (from MANE)
convertedCDNA <- convertedCDNA %>% rename(genomic_coordinateStart = start)
convertedCDNA <- convertedCDNA %>% rename(genomic_coordinateEnd = end)

# Set exon_number to start at 1 (important to calculate EJC position later)
convertedCDNA <- convertedCDNA %>% group_by(transcript_id) %>%
  mutate(exon_number = row_number()) %>%
  ungroup()

# ---------------------------------
# Deal with transcripts that do not map to Ensembl
transcript_subset <- convertedCDNA$transcript_id[ which(is.na(convertedCDNA$strand))]

# Calculate cDNA positions
threeUTR_intronSubset <- threeUTR_intronStop[ threeUTR_intronStop$transcript_id %in% transcript_subset, ]

threeUTR_intronSubset <- threeUTR_intronSubset %>% rename(genomic_coordinateStart = start)
threeUTR_intronSubset <- threeUTR_intronSubset %>% rename(genomic_coordinateEnd = end)

# Set exon_number to start at 1
threeUTR_intronSubset <- threeUTR_intronSubset %>% group_by(transcript_id) %>%
  mutate(exon_number = row_number()) %>%
  ungroup()

# Calculate width (end - start + 1)
threeUTR_intronSubset <- threeUTR_intronSubset %>%
  mutate(width = genomic_coordinateEnd - genomic_coordinateStart + 1)

# Calculate cumulative width
threeUTR_intronSubset <- as.data.table(threeUTR_intronSubset)
threeUTR_intronSubset[, cumwidth := cumsum(width), 
                      by = .(gene_name, transcript_id, type, strand)]

# Set the order by exon number
threeUTR_intronSubset[, exon_number := as.numeric(exon_number)]
setorder(threeUTR_intronSubset, transcript_id, exon_number)

# Update the exonStart_cDNA and exonEnd_cDNA positions
threeUTR_intronSubset[, exonStart_cDNA := cumwidth - width + 1]
threeUTR_intronSubset[, exonEnd_cDNA := cumwidth]

# Convert all values in "type" column to lowercase
# MANE uses "three_prime_UTR", Ensembl uses "three_prime_utr" - needs to be consistent
threeUTR_intronSubset$type <- tolower(threeUTR_intronSubset$type)

# Update the exonStart_cDNA for three_prime_utr rows
# Need exonStart_cDNA positions to begin from exonEnd_cDNA of stop_codon
for (i in 2:nrow(threeUTR_intronSubset)) {
  message(sprintf("%i of %i", i, nrow(threeUTR_intronSubset)))
  if(threeUTR_intronSubset$type[i] == "three_prime_utr"){
    prev_exon_end <- threeUTR_intronSubset$exonEnd_cDNA[i - 1]
    threeUTR_intronSubset$exonStart_cDNA[i] <- prev_exon_end + 1
    threeUTR_intronSubset$exonEnd_cDNA[i]   <- prev_exon_end + threeUTR_intronSubset$width[i]
  }
}

# Join with convertedCDNA data frame
setDT(convertedCDNA)
selected_columns <- c("transcript_id", "genomic_coordinateStart", "genomic_coordinateEnd", "exon_number",
                      "gene_name", "type", "strand", "cumwidth", "width", "exonStart_cDNA", "exonEnd_cDNA")

# Iterate through each selected column and update values in convertedCDNA
for(col in selected_columns){
  convertedCDNA[transcript_id %in% threeUTR_intronSubset$transcript_id,
                (col) := threeUTR_intronSubset[transcript_id %in% transcript_id, get(col)]]
}

# ---------------------------------
# Find introns between stop_codon and first 3' UTR exon

# Unique transcripts to iterate through
unique_transcripts <- unique(convertedCDNA$transcript_id)

# Iterate through each transcript
for(each in unique_transcripts){
  current_transcript <- convertedCDNA[transcript_id == each]
  
  if(current_transcript$strand[1] == "+"){
    if(nrow(current_transcript) > 1 &&
        current_transcript$genomic_coordinateStart[2] != current_transcript$genomic_coordinateEnd[1] + 1){
      
      new_row <- data.frame(
        transcript_id = each,
        genomic_coordinateStart = current_transcript$genomic_coordinateEnd[1] + 1,
        genomic_coordinateEnd = current_transcript$genomic_coordinateStart[2] - 1,
        gene_name = current_transcript$gene_name[1],
        type = "stop_intron",
        seqnames = current_transcript$seqnames[1],
        strand = current_transcript$strand[1],
        cumwidth = current_transcript$cumwidth[1],
        # width = start - end + 1, but using the adjacent start and ends (not the actual positions of the stop_intron)
        # Need to subtract 1  to get the actual width
        width = current_transcript$genomic_coordinateStart[2] - current_transcript$genomic_coordinateEnd[1] - 1,
        exon_number = 2
      )
      
      # Add the new stop_intron row between the stop_codon and first 3' UTR exon
      stop_codon_index <- which(convertedCDNA$transcript_id == each & convertedCDNA$type == "stop_codon")
      insert_position <- as.numeric(stop_codon_index)
      
      convertedCDNA <- rbind(convertedCDNA[1:insert_position, ], 
                             new_row, convertedCDNA[- (1:insert_position), ], fill = TRUE)
    }
  } # end of if loop for + strand
  
  else if(current_transcript$strand[1] == "-"){
    if(nrow(current_transcript) > 1 &&
       current_transcript$genomic_coordinateEnd[2] + 1 != current_transcript$genomic_coordinateStart[1]){
      
      new_row <- data.frame(
        transcript_id = each,
        genomic_coordinateStart = current_transcript$genomic_coordinateEnd[2] + 1,
        genomic_coordinateEnd = current_transcript$genomic_coordinateStart[1] - 1,
        gene_name = current_transcript$gene_name[1],
        type = "stop_intron",
        seqnames = current_transcript$seqnames[1],
        strand = current_transcript$strand[1],
        cumwidth = current_transcript$cumwidth[1],
        width = current_transcript$genomic_coordinateStart[1] - current_transcript$genomic_coordinateEnd[2] - 1,
        exon_number = 2
      )
      
      # Add the new stop_intron row between the stop_codon and first 3' UTR exon
      stop_codon_index <- which(convertedCDNA$transcript_id == each & convertedCDNA$type == "stop_codon")
      insert_position <- as.numeric(stop_codon_index)
      
      convertedCDNA <- rbind(convertedCDNA[1:insert_position, ], 
                             new_row, convertedCDNA[- (1:insert_position), ], fill = TRUE)
    }
  } # end of if loop for - strand
}

# Recalculate exon_number
convertedCDNA <- convertedCDNA %>% group_by(transcript_id) %>%
  mutate(exon_number = row_number()) %>%
  ungroup()

# Update the exonStart_cDNA and exonEnd_cDNA positions
setDT(convertedCDNA)
convertedCDNA[, cumwidth := cumsum(width), by = transcript_id]
convertedCDNA[, exonStart_cDNA := cumwidth - width + 1]
convertedCDNA[, exonEnd_cDNA := cumwidth]

# -------------------------------------------------------------------

write.table(convertedCDNA, file = sprintf("%s/threeUTRcDNA_intronStop.tsv", outDir),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)