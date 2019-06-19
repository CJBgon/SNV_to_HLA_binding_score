# Christian J Bouwens
# 2019-06-18
# visualise SNV - netMHCpan-4.0 results.

library(RColorBrewer)
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
countsum[,"ID"] <- factor(countsum$ID, levels = c("8mer", "9mer", "10mer", "peptidesum"))

meltdt <- melt(countsum, id = "ID")
display.brewer.all(colorblindFriendly = TRUE)
my_palette = brewer.pal(5, "Set2")[c(3,5,4)]
colour1 <- c("#67a9cf",
  "#999999",
  "#ffffff")

names(my_palette) <- levels(countsum$ID)[1:3]
colScale <- scale_colour_manual(name = "Var2",
                      values = my_palette,
                      aesthetics = "fill")


barresults <- ggplot(data = meltdt[!ID == "peptidesum"], aes(x = variable, y = value, fill = ID)) +
  geom_bar(stat = "identity")+
  labs(fill = "Peptide", x = "MHC-I", y = "Count")+
  scale_fill_brewer(palette = "Set2")

tiff(file = )
