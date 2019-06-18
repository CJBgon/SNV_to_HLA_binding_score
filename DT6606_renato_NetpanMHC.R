#!/usr/bin/Rscript

# Christian J Bouwens
# 2019-06-13
# Script that returns 8, 9 and 10mer peptides containing SNV from
# a .vcf or .tsv file of single nucleotide variations from whole exome
# sequencing analysis.

suppressPackageStartupMessages(library(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
library(UniProt.ws)
library(data.table)
library(biomaRt)


option_list = list(
  make_option(c("-f", "--file"),
   action="store",
   default=NA,
   type='character',
   help="location of the .vcf or .tsv file"),
  make_option(c("-l", "--length"),
   action="store",
   default=NA,
   type='numeric',
   help="length (rows) of the VCF header.")
)
opt = parse_args(OptionParser(option_list = option_list))


#functions
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


seqchecking <- function(gene, loc, from, seqdata,
  namecol = "V1", seqcol = "V2") {
  # Checks if the a gene name and expected amino acid are in one or more
  # of the sequences provided.
  #
  # Args:
  #   gene: Gene name of interest. e.g.: "KRAS"
  #   loc: location of the expected amino acid. e.g.: 12 in G12D.
  #   from: expected amino acid in sequence, upper case single letter annot.
  #         "G" in G12D.
  #   seq: Data frame with names and sequences that should be checked.
  #   namecol: name of the column containing gene names in seq. (default = "V1")
  #   seqcol: name of the column containing the peptide sequences in
  #         seq. (default = "V2")
  #
  # Returns:
  #   A vector with sequences that fit the gene name and expected AA.
  library(data.table)

  if (!is.data.table(seqdata)){
   seqdataDT <- as.data.table(seqdata)
  }

  if (is.na(gene) || is.null(gene)){
    stop("please enter a valid gene name.")
  } else if (gene %in% unlist(seqdataDT[,get(namecol)])) {
    sequence <-
      unlist(seqdataDT[get(namecol) == gene, get(seqcol)], use.names = F)
  } else{sequence <- NA}

  returnseq <- c()
  for (i in sequence) {
    if (is.na(i)) {
      returnseq <- c(returnseq, NA)
    }else if(substring(i, first = loc, last = loc) == from) {
      returnseq <- c(returnseq, i)
    }else{returnseq <- c(returnseq, NA)}
  }
  return(returnseq)
}


seqselect <- function(seqvector) {
 # filters NA elements and selects the first element out of a
 # list of multiple choices per element i.
 #
 # Args:
 #    seqvector: an list containing NA and sequences.
 #
 #Returns: A vector of single elements extracted from seqvector. Returns NA if
 #    there was no element to extract or if the only option was NA.
 #
 #
 #
 select <-c()
  for (i in seq_along(seqvector)) {
    narem <- c()
    if (length(seqvector[[i]]) == 1){
      select <- c(select, seqvector[[i]])

    }else if (length(seqvector[[i]]) >= 2){

      narem <- seqvector[[i]][!is.na(seqvector[[i]])]

      if(is.null(narem) || length(narem) <1) {
        select <- c(select, NA)
      }else{ select <- c(select, narem[1])}
    }else{select <- c(select, NA)}
  }
  return(select)
}


formatmassager <- function (data, out, sep = "", remove1 = NULL,
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
  #
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

# remove .VCF annotation data in the top rows,
# the length of these might differ depending on the pipeline used.
# to do: give an optoin to either work with VCF or TSV/CSV files.
usedat <- dat[ -(1:opt$l)]  #opt$l
# gather nonsynonymous, exonic variations.
nonsynonmut <- usedat[ExonicFunc.refGene == "nonsynonymous SNV"]

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

# additional biomart query
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")

# because the UNIprot ID does no always net a sequence we go over the mutations
# again but query the refseq ID and otherwise the protein name.
# which makes it so unless we query individually,
# we cant see the differences between genen names
# for now we will use the funciton nm_seqchecked to get sequences
# that have the expected aminoacid and gene name to match the querried
# peptide sequence to the mutation in the VCF/TSV file.

 # 1. remove get dataframe_NA and dataframe 1.
naAAdata <- AAdata[unlist(is.na(AAdata[, 11]), use.names = FALSE),]
valAAdata <- AAdata[!unlist(is.na(AAdata[, 11]), use.names = FALSE),]

 # 2. query dataframe_NA for refeseq peptide.
martseq <- getSequence(id = unlist(naAAdata[, 2], use.names = F),
   type = "refseq_mrna",
   seqType = "peptide",
   mart = mart,
   verbose = F)

 # 3. check which refseq hits match using seqchecking.
 # append those to dataframe 1.
nm_seqchecked <-
 mapply(seqchecking,
   gene = naAAdata[, 2],
   loc = naAAdata[, 6],
   from = naAAdata[, 4],
   MoreArgs = list(
   seqdata = martseq,
   namecol = "refseq_mrna",
   seqcol = "peptide")
 )

nm_seq_selected <- seqselect(nm_seqchecked)
for (i in seq_along(nm_seq_selected)) {
  naAAdata[i,11] <- nm_seq_selected[i]
}

 # 4. for the hits that dont match (including still NA) check on mgi_symbol.
naAAdata2 <- naAAdata[unlist(is.na(naAAdata[, 11]), use.names = F),]
val_naAAdata <- naAAdata[!unlist(is.na(naAAdata[, 11]), use.names = F),]

martseq_mgi <- getSequence(id = unlist(naAAdata[, 1], use.names = F),
   type = "mgi_symbol",
   seqType = "peptide",
   mart = mart,
   verbose = F)

mgi_seqchecked <-
mapply(seqchecking,
  gene = naAAdata2[, 1],
  loc = naAAdata2[, 6],
  from = naAAdata2[, 4],
  MoreArgs = list(
  seqdata = martseq_mgi,
  namecol = "mgi_symbol",
  seqcol = "peptide")
)

mgi_seq_selected <- seqselect(mgi_seqchecked)
for (i in seq_along(mgi_seq_selected)) {
naAAdata2[i,11] <- mgi_seq_selected[i]
}

 # 5. append after seqchecking and place the remaining NAs in a seperate file.
naAAdata2_2 <-naAAdata2[unlist(is.na(naAAdata2[, 11]), use.names = F),]
val_naAAdata2 <- naAAdata2[!unlist(is.na(naAAdata2[, 11]), use.names = F),]


AAdata_fin<- rbind(valAAdata, val_naAAdata, val_naAAdata2)
naAAdata_excluded <- naAAdata2_2

write.csv(naAAdata_excluded,
  file = "/home/max/Documents/GitHub/gits/excluded_mut.csv")

eightmers <- mapply(Ncharextract,
  AAdata_fin[, 11],
  n = 8,
  loc = AAdata_fin[, 6],
  from = AAdata_fin[, 4],
  mut = AAdata_fin[, 5],
  replace = TRUE
  )

ninemers <- mapply(Ncharextract,
  AAdata_fin[, 11],
  n = 9,
  loc = AAdata_fin[, 6],
  from = AAdata_fin[, 4],
  mut = AAdata_fin[, 5],
  replace = TRUE
  )

tenmers <- mapply(Ncharextract,
  AAdata_fin[, 11],
  n = 10,
  loc = AAdata_fin[, 6],
  from = AAdata_fin[, 4],
  mut = AAdata_fin[, 5],
  replace = TRUE
  )

AAdata_fin2 <- cbind(AAdata_fin, eightmers, ninemers, tenmers)
write.csv(AAdata_fin2,
file = "/home/max/Documents/GitHub/gits/AA_peptides_mutations.csv")

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
  file = "/home/max/Documents/GitHub/gits/DT6606_eightmer.txt",
  quote = FALSE,
  sep = "",
  row.names = FALSE,
  col.names = FALSE)

write.table(x = unlist(ninemerstrings),
  file = "/home/max/Documents/GitHub/gits/DT6606_ninemer.txt",
  quote = FALSE,
  sep = "",
  row.names = FALSE,
  col.names = FALSE)

write.table(x = unlist(tenmerstrings),
  file = "/home/max/Documents/GitHub/gits/DT6606_tenmer.txt",
  quote = FALSE,
  sep = "",
  row.names = FALSE,
  col.names = FALSE)
