#~!/bin/bash
export PATH=~/Documents/netMHC-4.0:$PATH
export PATH=$PATH:/usr/bin/Rscript
# run the R script DT6606_renato_NetpanMHC.R with the VCF file as variable.
chmod +x /home/christian/Documents/GitHub/gits/DT6606_renato_NetpanMHC.R

cd /home/christian/Documents/GitHub/gits
./DT6606_renato_NetpanMHC.R -f /home/christian/Documents/GitHub/gits/DT6606_variants_sub.tsv -l 30

cd ../
#run the netpanMHC script

netMHC -a H-2-Db -f /home/christian/Documents/GitHub/gits/DT6606_eightmer.txt -p -l 8 -xls -xlsfile /home/christian/Documents/GitHub/gits/NETPAN_Db_eightmer.xls
netMHC -a H-2-Kb -f /home/christian/Documents/GitHub/gits/DT6606_eightmer.txt -p -l 8 -xls -xlsfile /home/christian/Documents/GitHub/gits/NETPAN_Kb_eightmer.xls

netMHC -a H-2-Db -f /home/christian/Documents/GitHub/gits/DT6606_ninemer.txt -p -l 9 -xls -xlsfile /home/christian/Documents/GitHub/gits/NETPAN_Db_ninemer.xls
netMHC -a H-2-Kb -f /home/christian/Documents/GitHub/gits/DT6606_ninemer.txt -p -l 9 -xls -xlsfile /home/christian/Documents/GitHub/gits/NETPAN_Kb_ninemer.xls

netMHC -a H-2-Db -f /home/christian/Documents/GitHub/gits/DT6606_tenmer.txt -p -l 10 -xls -xlsfile /home/christian/Documents/GitHub/gits/NETPAN_Db_tenmer.xls
netMHC -a H-2-Kb -f /home/christian/Documents/GitHub/gits/DT6606_tenmer.txt -p -l 10 -xls -xlsfile /home/christian/Documents/GitHub/gits/NETPAN_Kb_tenmer.xls
