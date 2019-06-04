if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UniProt.ws")
BiocManager::install("biomaRt")
install.packages("data.table")
install.packages("styler")


library(UniProt.ws)
library(data.table)
library(biomaRt)

# load in annotated whole exsome sequencing single nucleotide variation data.
read.table("~/Documents/GitHub/gits/DT6606_variants_sub.tsv",
  sep = "\t",
  fill = T)

dat <- fread(file = "~/Documents/GitHub/gits/DT6606_variants_sub.tsv",
  sep2 = ";",
  fill = F,
  na.strings ="NA",
  check.names = T,
  header = T)
# remove .VCF data in the first top 30 rows.
# the length of these might differ depending on the pipeline used.
# to do: give an optoin to either work with VCF or TSV/CSV files.
usedat <- dat[ -(1:30)]
# gather nonsynonymous, exonic variations.
nonsynonmut <- usedat[ExonicFunc.refGene == "nonsynonymous SNV"]



rsubstr <- function(x,n) {
# substring only works from left to right.
# this function grabs n characters from the end of the string.

  substring(x, nchar(x) - n + 1)
}

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

# TO DO: add maxlength to give the maximum number of cycles per mutation.
#        create while (sequence is NA AND k> max number of cycles, k +5 (to the new NM_####)) loop and run again.
#lapply(aasplit, function(x) {
#  unit1 <- unlist(x)[2] #hier neemt hij alleen de eerste van elke list.
#  k=2
#  sapply(unit1, function(i) {
  #  match(i[2], )
  #}
#})


#firstele <- lapply(aasplit, "[[", 1)

#aalistlengths <- lapply(AAmutsplit, length)

#AAmutsplit[[8]][2]
#class(AAmutsplit)
#str(AAmutsplit)
#length(AAmutsplit)
#AAmustsep <- sub(x = AAmutsmall,
#  replacement = "",
#  pattern = ".*,")sapply(AAmutsplit, function(x) {



# grab the NCBI refseq reference, gene name, original amino acid,
# location of the mutation and mutation.
NM  <- unlist(dfaasplit[, 3], use.names = F)
mut <- unlist(dfaasplit[, 6], use.names = F)
mutpos <- as.numeric(gsub(".*?([0-9]+).*",
  "\\1",
  x = mut))
from <- substr(mut, 3, 3)
to <- sapply(mut, function(x) {
  substring(x, nchar(x))
  }
  , USE.NAMES = F)
name <- unlist(dfaasplit[1, ], use.names = F)

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

NMmatch <- match(unlist(dfaasplit[2 ,]),resref_fin[, 1],)
AAmat <- resref_fin[NMmatch, c(-3, -4)]
AAdata <- cbind(name,
  NM,
  mut,
  from,
  to,
  mutpos,
  AAmat)
#end-----

# biomart approach
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
listAttributes(mart)
searchAttributes(mart, pattern = "symbol")
searchAttributes(mart, pattern = c("peptide"))
searchFilters(mart, pattern = "mgi")
annot <- getBMlist(attributes =
  "ensembl_gene_id",
  "refseq_mrna",
  "peptide",
  filters = "refseq_mrna",
  values  = NM,
   mart=mart)

martseq_nm <- getSequence(id = NM,
  type = "refseq_mrna",
  seqType = "peptide",
  mart = mart,
  verbose = FALSE)
# how many will return if I use simple gene names?
martseq_mgi <- getSequence(id = name,
  type = "mgi_symbol",
  seqType = "peptide",
  mart = mart)

length(unique())
length(unique(NM))#30 unique refseq codes
length(unique(martseq_nm[, 2])) #24/30
length(unique(martseq_mgi[, 2])) #17/30 met mgi_symbol
length(unique(resref_fin[,7])) #19/30 met uniprotID.

# using names doesn't really work because of transcript variants.
# we should annotate into ENSGM and ENSTM. those should always net the correct
# protein sequence.
martseqmgi2 <- martseq_mgi[!martseq_mgi[, 1] == "Sequence unavailable", ]

#lets continue with the resref_fin sequences.

resref_fin


Ncharextract <- function(x, n, loc, mut = "X") {
  # extracts +n and -n characters surrouding a location in a vector.
  #
  # Args:
  #   x: A vector where the characters will be extracted from
  #   n: The amount of characters to be extracted from the vect,
  #      should be a positive integer.
  #   loc: Numeric location in the vector from where
  #      surrounding characters are extracted.
  #   replace: If TRUE, replaces the character at position loc
  #      with the character supplied in mut. Default is TRUE
  #   mut: character that replaces original character in location loc.
  #
  # Returns:
  #   A vector of 2n+1 length extracted from vector x.
  mer <- c()
  if (is.na(x)) {
    stop("x is not a character or numeric vector, or x is not supplied.")
  }else{
     Rmer <- substring(x, first = loc - n, last = loc - 1)
     Lmer <- substring(x, first = loc + 1, last = loc + n)
  }

  if (replace == TRUE) {
    mid <- mut
  }else{
    mid <- substring(x, first = loc, last = loc)
  }

  mer <- c(Rmer + mid + Lmer)
}


UniproteinQC <- function (data, data2, seqcol, loc, expectedAA){
  # qualit control function that checks if the the protein sequence was
  # successfully found and if the expected amino acid is in the right place.
  #
  # Args:
  #   data: Data.table with protein names, extracted refseq IDs, Uniprot IDs
  #         and protein sequences.
  #   data2: Initial VCF or .tsv data containing the CNV and annotation.
  #   seqcol: column in data containing the protein sequences.
  #   loc: column in data containing the mutation site.
  #
  # Returns:
  #   A revised data.table with adjusted protein sequences where previously
  #   sequences were missing or incorrect.
  for (i in seq(along = data)) {

    if (is.na(AAdata[i,seqcol])){ #  first attempt to resolve all NA sequences.

    }
    # subsequently check if found protein sequences contain the expected AA.
}


#check for each protein sequence if the expected AA(from) is in position mut.
mer8 <- c()
mer9 <- c()
mer10 <- c()

for(i in seq_along(AAdata[,1])){

  if (is.na(AAdata[i,"SEQUENCE"])){
    mer8 <- c(mer8, "nosequence")
    mer9 <- c(mer9, "nosequence")
    mer10 <- c(mer10, "nosequence")
  }else if (substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"], last = AAdata[i,"mutpos"]) == AAdata[i,"from"]){
    #select -7 AA
    rightmer8 <- c(substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"]-7, last = AAdata[i,"mutpos"]-1))
    #select -8AA
    rightmer9 <- c(substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"]-8, last = AAdata[i,"mutpos"]-1))
    #select -9AA
    rightmer10 <- c(substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"]-9, last = AAdata[i,"mutpos"]-1))
    #select +7 AA
    leftmer8 <- c(substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"]+1, last = AAdata[i,"mutpos"]+7))
    #select +8AA
    leftmer9 <- c(substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"]+1, last = AAdata[i,"mutpos"]+8))
    #select +9AA
    leftmer10 <- c(substring(AAdata[i,"SEQUENCE"], first = AAdata[i,"mutpos"]+1, last = AAdata[i,"mutpos"]+9))

    #stitch them together
    comb8 <- paste(rightmer8, AAdata[i,"to"], leftmer8, sep = "")
    comb9 <- paste(rightmer9, AAdata[i,"to"], leftmer9, sep = "")
    comb10 <- paste(rightmer10, AAdata[i,"to"], leftmer10, sep = "")

    mer8 <- c(mer8, comb8)
    mer9 <- c(mer9, comb9)
    mer10 <- c(mer10, comb10)

  }else{
    mer8 <- c(mer8, "AA not found in sequence")
    mer9 <- c(mer9, "AA not found in sequence")
    mer10 <- c(mer10, "AA not found in sequence")
  }

}
AAdata2 <- cbind(AAdata, mer8, mer9, mer10)

#WIP find a way to deal with transcript variants when they don't have an uniprot sequence attached.


write.csv(x=AAdata2,file = "~/Documents/MHC_neoantigen_select_WES/DT6606_renato_NetpanMHC.csv", row.names = F, col.names = T)


#filter rout unreviewed sequences:
dat <- fread("~/Documents/MHC_neoantigen_select_WES/DT6606_renato_NetpanMHC.csv")
rev_dat <- dat[REVIEWED == "reviewed" & !is.na(SEQUENCE)]
write.csv(x=AAdata2, file= "~/Documents/MHC_neoantigen_select_WES/DT6606_renato_NetpanMHC_filtered.csv")
