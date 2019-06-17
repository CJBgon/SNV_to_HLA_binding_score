#!/usr/bin/Rscript

# Christian J Bouwens
# 2019-06-13
# Script that returns 8, 9 and 10mer peptides containing SNV from
# a .vcf or .tsv file of single nucleotide variations from whole exsome
# sequencing analysis.

suppressPackageStartupMessages(library(optparse)) # don't say "Loading required package: optparse"
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf
library(UniProt.ws)
library(data.table)
library(biomaRt)



option_list = list(
  make_option(c("-f", "--file"), action="store", default=NA, type='character',
              help="location of the .vcf or .tsv file"),
  make_option(c("-l", "--length"), action="store", default=NA, type='numeric',
              help="length of of the VCF header. ")
)
opt = parse_args(OptionParser(option_list=option_list))

rsubstr <- function(x,n) {
  # substring only works from left to right.
  # this function grabs n characters from the end of the string.

  substring(x, nchar(x) - n + 1)
}


stringslider <- function(stri, l) {
  # Returns a vector of unique strings of length l extracted from string
  # works through a sliding window from left to right. Checks if each returned
  # string is of length l.
  #
  # Args:
  #   stri: A string from on which the extraction occurs.
  #   l: length of the extracted strings.
  #
  # Returns:
  #   A vector of extracted strings with length l.
  strle <- nchar(stri)
  i <- 1
  subc <- c()

  while (i + l - 1 <= strle) {
    subc <- c(subc,
      substr(stri,
       start = i,
       stop  = i + l - 1))
   i <- i + 1
  }
  return(subc)
}


Ncharextract <- function(x, n, loc, from, replace = TRUE, mut = "X") {
  # extracts +n and -n characters surrouding a location in a vector.
  #
  # Args:
  #   x: A vector where the characters will be extracted from
  #   n: The amount of characters to be extracted from the vect,
  #      should be a positive integer.
  #   loc: Numeric location in the vector from where
  #      surrounding characters are extracted.
  #   from : Original amino acid in the mutated location. e.g. A in A12K.
  #   replace: If TRUE, replaces the character at position loc
  #      with the character supplied in mut. Default is TRUE
  #   mut: character that replaces original character in location loc.
  #
  # Returns:
  #   A vector of 2n+1 length extracted from vector x.
  mer <- c()
  if (is.na(x)) {
    stop("x is not a character or numeric vector, or x is not supplied.")
  }else if (x == "Sequence unavailable") {
    mer <- NA
  }else if (0.5 * n > nchar(x)) {
    mer <- "Sequence too short"
  }else if (substring(x, first = loc, last = loc) == from) {
     Rmer <- substring(x, first = loc - n, last = loc - 1)
     Lmer <- substring(x, first = loc + 1, last = loc + n)

  if (replace == TRUE) {
    mid <- mut
  }else{
    mid <- substring(x, first = loc, last = loc)
  }

  mer <- paste(Rmer, mid, Lmer, sep = "")
  }else{
  mer <- "Incorrect mutation location or sequence"
  }
}


formatmassager <- function (data, out, sep ="", remove1 = NULL,
  remove2 = NULL, savetxt = TRUE){
  # A vector of short peptide sequences (8, 9 and 10 mers
  # and format them for netpanMHC analysis.
  #
  # Args:
  #   data: Data.table or vector file with character strings.
  #   out: location and filename for the .peptide file.
  #   remove1: "character" that should be filtered out of the vector.
  #   remove2: "character" that should be filtered out of the vector.
  #   savetxt: save the data.table as a .txt file? logical, default = TRUE
  # Returns:
  #   saves a .txt file with newline delimited of unique character strings from
  #   the vector x.
  res <- unique(data[!is.na(data) & !data %in% c(remove1, remove2)])
  if (savetxt == TRUE) {
  write.table(x = res,
    file = out,
    quote = FALSE,
    sep = sep,
    row.names = FALSE,
    col.names = FALSE)
  }
  return(res)
}


dat <- fread(file = opt$f,
  sep = "\t",
  sep2 = ";",
  fill = T,
  na.strings ="NA",
  check.names = T,
  header = T)

# remove .VCF data in the first top 30 rows.
# the length of these might differ depending on the pipeline used.
# to do: give an optoin to either work with VCF or TSV/CSV files.
usedat <- dat[ -(1:opt$l)]
# gather nonsynonymous, exonic variations.
nonsynonmut <- usedat[ExonicFunc.refGene == "nonsynonymous SNV"]

# multiple NCBI gene variants and corresponding CNV location
# are located in the same column, we extract the last 50 digits and stringsplit
# this the final reference and mutation are available.
# TO DO: grab all combinations and use those downstream in case no sequence and
# AA can be found for the first reference in the UniProt database.
AAmutsmall <- nonsynonmut[, sub(x = AAChange.refGene,
  replacement = "",
  pattern = ".*(?=.{50}$)",
  perl = T)]

# We take all NM options and locations per SNV
# While this might increase computational cost we don't need to run a loop to
# check if every refseq ID has resulted in a protein sequence
# This increases the change of finding the correct protein sequence in UniProt
# At the end we'll filter for unique peptides
AAmutsplit <- strsplit(nonsynonmut[, AAChange.refGene],
  ",",
  perl = T)

# The AAChange.refGene are separated by ":", we need each in a seperate column
aasplit <- lapply(AAmutsplit, function(x) {
  strsplit(x, split = ":")
})

# turn into a data.frame instead of data.table or unlist as we need each list
# element as a seperate column
dfaasplit <- as.data.frame(aasplit,
  stringsAsFactors = F)

dfaasplit <- t(dfaasplit)

# grab the NCBI refseq reference, gene name, original amino acid,
# location of the mutation and mutation.
NM  <- unlist(dfaasplit[, 2], use.names = F)
mut <- unlist(dfaasplit[, 5], use.names = F)
mutpos <- as.numeric(gsub(".*?([0-9]+).*",
  "\\1",
  x = mut))
from <- substr(mut, 3, 3)
to <- sapply(mut, function(x) {
  substring(x, nchar(x))
  }
  , USE.NAMES = F)
name <- unlist(dfaasplit[, 1], use.names = F)

# select m. musculus uniprot datbase for AA sequence extraction.----
mouseUp <- UniProt.ws(10090)
# determine all refseq names
mouseke <- keys(mouseUp, "REFSEQ_NUCLEOTIDE")
# remove .# version info
mouseke2 <- sub(x = mouseke,
   pattern = "\\..*",
   replacement = "")
# match our NM_######## linked to non-synonymous mutations to the refseq keys
keymou <- match(NM, mouseke2)
# select keys of each mutation as they appear in the uniprot dataset.
# this step is neccesary as we are missing the version info.
keymouref <- mouseke[keymou]
columnref <- c("UNIPROTKB")

# get uniprotKB identifier
# as not all refseq identfiers are linked to a sequence
# some refseq nucleotides are converted to different UNIPROTKB identifers
# in this step we idetify these cases
res <- select(x = mouseUp,
  keys = keymouref,
  columnref,
  keytype = "REFSEQ_NUCLEOTIDE")

# we can then extract the correct UNIPROT id
columns <- c( "REVIEWED", "SEQUENCE")
uniprot <- res[,"UNIPROTKB"]
protsub <- sub(uniprot,
  pattern = "(.*->)",
  replacement = "",
  perl =T)

protsubsub <- sub(protsub,
  pattern = "(-.*)",
  replacement = "",
  perl = T)

resref <- cbind(res, protsub)

# select sequences of eacch uniprot IDs.
res2 <- select(x=mouseUp,
  keys = resref[, 3],
  columns,
  keytype = "UNIPROTKB")

# trace back uniprotID linked sequences with our NM with.
NMsub <- sub(x = resref[, 1],
  pattern = "\\..*",
  replacement = "")

# build up final table
resref_fin <- cbind(NMsub, resref, res2)

NMmatch <- match(unlist(dfaasplit[, 2],
  use.names = F),
  resref_fin[, 1])

AAmat <- resref_fin[NMmatch, c(-3, -4)]
AAdata <- cbind(name,
  NM,
  mut,
  from,
  to,
  mutpos,
  AAmat,
  stringsAsFactors = F)
# end of uniprot query----
AAdata
# additional biomart query
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")

# because the UNIprot ID does no always net a sequence we go over the mutations
# again but query the refseq ID and otherwise the protein name.
for (i in seq_along(AAdata[,1])) {
  # If there is no sequence (NA), search biomaRt for the refeseq sequence.
  if (is.na(AAdata[i, 11])) {
    martseq_nm <- getSequence(id = unlist(AAdata[i, 2], use.names = F),
      type = "refseq_mrna",
      seqType = "peptide",
      mart = mart,
      verbose = F)

      # If there is no sequence linked to the refseq id, try the  gene name.
      if (is.na(martseq_nm[1, 1]) ||
       martseq_nm == "Sequence unavailable"){
        martseq_mgi <- getSequence(id = unlist(AAdata[i, 1], use.names = F),
          type = "mgi_symbol",
          seqType = "peptide",
          mart = mart)

          # if martseq mgi does not return NA or Sequence unavailable or is not null use it.
          if(!is.na(martseq_mgi[1, 1]) ||
          martseq_mgi[1, 1] != "Sequence unavailable" ||
          !is.null(martseq_mgi)) {
            AAdata[i, 11] <- unlist(martseq_mgi[1, 1], use.names = F)
          }
      }else{
        AAdata[i, 11] <- unlist(martseq_nm[1, 1], use.names = F)
      }
    }
  }


eightmers <- mapply(Ncharextract,
  AAdata[, 11],
  n = 8,
  loc = AAdata[, 6],
  from = AAdata[, 4],
  mut = AAdata[, 5],
  replace = TRUE
  )

ninemers <- mapply(Ncharextract,
  AAdata[, 11],
  n = 9,
  loc = AAdata[, 6],
  from = AAdata[, 4],
  mut = AAdata[, 5],
  replace = TRUE
  )

tenmers <- mapply(Ncharextract,
  AAdata[, 11],
  n = 10,
  loc = AAdata[, 6],
  from = AAdata[, 4],
  mut = AAdata[, 5],
  replace = TRUE
  )

AAdata2 <- cbind(AAdata, eightmers, ninemers, tenmers)


# create the textfiles for netpanMHC scoring:
eightmersclean <- formatmassager(data = unlist(eightmers),
  savetxt = FALSE,
  remove1 = "Incorrect mutation location or sequence",
  remove2 = "Sequence too short")

ninemersclean <- formatmassager(data = unlist(ninemers),
  savetxt = FALSE,
  remove1 = "Incorrect mutation location or sequence",
  remove2 = "Sequence too short")

tenmersclean <- formatmassager(data = unlist(tenmers),
  savetxt = FALSE,
  remove1 = "Incorrect mutation location or sequence",
  remove2 = "Sequence too short")

eightmerstrings <- mapply(stringslider, eightmersclean, 8, USE.NAMES=FALSE)
ninemerstrings  <- mapply(stringslider, ninemersclean,  9, USE.NAMES=FALSE)
tenmerstrings   <- mapply(stringslider, tenmersclean,   10, USE.NAMES=FALSE)

write.table(x = unlist(eightmerstrings),
  file = "/home/christian/Documents/GitHub/gits/DT6606_eightmer.txt",
  quote = FALSE,
  sep = "",
  row.names = FALSE,
  col.names = FALSE)

write.table(x = unlist(ninemerstrings),
  file = "/home/christian/Documents/GitHub/gits/DT6606_ninemer.txt",
  quote = FALSE,
  sep = "",
  row.names = FALSE,
  col.names = FALSE)

write.table(x = unlist(tenmerstrings),
  file = "/home/christian/Documents/GitHub/gits/DT6606_tenmer.txt",
  quote = FALSE,
  sep = "",
  row.names = FALSE,
  col.names = FALSE)
