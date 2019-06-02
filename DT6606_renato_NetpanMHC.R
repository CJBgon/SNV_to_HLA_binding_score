if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UniProt.ws")
install.packages("data.table")
install.packages("styler")
install.packages("lintr")
library(UniProt.ws)
library(data.table)
library(lintr)

# load in annotated whole exsome sequencing single nucleotide variation data.
read.table("~/Documents/GitHub/gits/DT6606_variants_sub.tsv", sep = "\t", fill = T)
dat <- fread(file = "~/Documents/GitHub/gits/DT6606_variants_sub.tsv",
  sep2 = ";",
  fill = F,
  na.strings ="NA",
  check.names = T)
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

#the
AAmutsmall <- nonsynonmut[ ,sub(x = AAChange.refGene,
  replacement = "",
  pattern = ".*(?=.{50}$)",
  perl = T)]
AAmustsep <- sub(x = AAmutsmall,
  replacement = "",
  pattern = ".*,")
AAsplit <- strsplit(x = AAmustsep, split = ":")
AAsplit <- as.data.table(AAsplit)

NM <- unlist(AAsplit[2, ], use.names = F)
mut <- unlist(AAsplit[5, ], use.names = F)
mutpos <- as.numeric(gsub(".*?([0-9]+).*", "\\1", x = mut))
from <- substr(mut, 3,3)
to <- sapply(mut, function(x){
  substring(x,nchar(x))
}, USE.NAMES = F)
name <- unlist(AAsplit[1,], use.names=F)

#select m. musculus uniprot datbase for AA sequence extraction.
mouseUp <- UniProt.ws(10090)
columns(mouseUp)
mouseke <- keys(mouseUp, "REFSEQ_NUCLEOTIDE")#determine all refseq names
mouseke2 <- sub(x = mouseke, pattern = "\\..*", replacement = "") #remove .# version info
keymou <- match(NM,mouseke2) #match our NM_######## linked to non-synonymous mutations to the reseq keys
keymouref <- mouseke[keymou] #select keys of each mutation

columnref <- c("UNIPROTKB")
res <- select(x = mouseUp, keys = keymouref , columnref, keytype =  "REFSEQ_NUCLEOTIDE")#get uniprotKB identifier (as not all refseq identfiers are linked to a sequence
columns <- c( "REVIEWED", "SEQUENCE")
uniprot <- res[,"UNIPROTKB"]
protsub <- sub(x=uniprot, pattern = "(.*->) ", replacement = "", perl =T)
protsubsub <- sub(x=protsub, patter = "(-.*)", replacement ="", perl = T)

resref <- cbind(res,protsubsub)
res2 <- select(x=mouseUp, keys = resref[,3], columns, keytype = "UNIPROTKB") #select sequences of eacch uniprot IDs.
#trace back uniprotID linked sequences with our NM with.
NMsub <- sub(x = resref[,1], pattern = "\\..*", replacement = "")
resref_fin <- cbind(NMsub, resref, res2)

NMmatch <- match(NM,resref_fin[,1])
AAmat <- resref_fin[NMmatch,c(-3, -4)]
AAdata <- cbind(name, NM, mut, from, to, mutpos, AAmat)

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
