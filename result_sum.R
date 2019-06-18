# Christian J Bouwens
# 2019-06-18
# visualise SNV - netMHCpan-4.0 results.

library(ggplot2)
library(data.table)

# load results.
eightdb <- fread("~/NGS/Blood_Renato_DT6606_WES/neoantigen_prediction/NETPAN_Db_eightmer.xls")
ninedb <- fread("~/NGS/Blood_Renato_DT6606_WES/neoantigen_prediction/NETPAN_Db_ninemer.xls")
tendb <- fread("~/NGS/Blood_Renato_DT6606_WES/neoantigen_prediction/NETPAN_Db_tenmer.xls")
eightkb <- fread("~/NGS/Blood_Renato_DT6606_WES/neoantigen_prediction/NETPAN_Kb_eightmer.xls")
ninekb <- fread("~/NGS/Blood_Renato_DT6606_WES/neoantigen_prediction/NETPAN_Kb_ninemer.xls")
tenkb <- fread("~/NGS/Blood_Renato_DT6606_WES/neoantigen_prediction/NETPAN_Kb_tenmer.xls")

listmers <- list(eightdb,
  ninedb,
  tendb,
  eightkb,
  ninekb,
  tenkb)
# Select low and high binders@
# cutoffs determined by netMHCpan are:
# low: 2.00 < x > 0.50 Rank Threshold
# high: x < 0.50 Rank Threshold

lowbinders <- lapply(listmers, function(x) {
  x[Rank >= 0.5 & Rank < 2, ]
})

highbinders <- lapply(listmers, function(x) {
  x[Rank <= 0.5]
})

# questions to answer: how many haplotype? -barplot
dbcounthigh <- lapply(highbinders[1:3], nrow)
kbcounthigh <- lapply(highbinders[4:6], nrow)
dbcountlow <- lapply(lowbinders[1:3], nrow)
kbcountlow  <- lapply(lowbinders[4:6], nrow)

countdt <- data.table(c("8mer", "9mer", "10mer"),
  unlist(dbcounthigh),
  unlist(kbcounthigh),
  unlist(dbcountlow),
  unlist(kbcountlow))

countdt[, sum:= V2 + V3 + V4 + V5]
rowcount <- countdt[, lapply(.SD,sum), .SDcols = colnames(countdt[,2:6])]
rowcount <- cbind(V1 = "peptidesum", rowcount)
countsum <- rbind(countdt, rowcount)
colnames(countsum) <- c("ID", "Db_High", "Kb_High", "Db_low", "Kb_low", "sum")
