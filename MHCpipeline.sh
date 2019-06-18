#~!/bin/bash
export PATH=/home/max/Documents/MHC_neoantigen_select_WES:$PATH
export PATH=$PATH:/usr/bin/Rscript
# run the R script DT6606_renato_NetpanMHC.R with the VCF file as variable.
chmod +x /home/max/Documents/GitHub/gits/variantstopeptides.R

cd /home/max/Documents/GitHub/gits
./variantstopeptides.R -f /home/max/NGS/Blood_Renato_DT6606_WES/variant_calling/blood_SNP_calling.Somatic.hc.mm10_multianno -l 30

cd ../
#run the netpanMHC script

cd /home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0

/home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0/netMHCpan -a H2-Db -f /home/max/Documents/GitHub/gits/DT6606_eightmer.txt -p -l 8 -xls -xlsfile /home/max/Documents/GitHub/gits/NETPAN_Db_eightmer.xls
/home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0/netMHCpan -a H2-Kb -f /home/max/Documents/GitHub/gits/DT6606_eightmer.txt -p -l 8 -xls -xlsfile /home/max/Documents/GitHub/gits/NETPAN_Kb_eightmer.xls

/home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0/netMHCpan -a H2-Db -f /home/max/Documents/GitHub/gits/DT6606_ninemer.txt -p -l 9 -xls -xlsfile /home/max/Documents/GitHub/gits/NETPAN_Db_ninemer.xls
/home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0/netMHCpan -a H2-Kb -f /home/max/Documents/GitHub/gits/DT6606_ninemer.txt -p -l 9 -xls -xlsfile /home/max/Documents/GitHub/gits/NETPAN_Kb_ninemer.xls

/home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0/netMHCpan -a H2-Db -f /home/max/Documents/GitHub/gits/DT6606_tenmer.txt -p -l 10 -xls -xlsfile /home/max/Documents/GitHub/gits/NETPAN_Db_tenmer.xls
/home/max/Documents/MHC_neoantigen_select_WES/netMHCpan-4.0/netMHCpan -a H2-Kb -f /home/max/Documents/GitHub/gits/DT6606_tenmer.txt -p -l 10 -xls -xlsfile /home/max/Documents/GitHub/gits/NETPAN_Kb_tenmer.xls
