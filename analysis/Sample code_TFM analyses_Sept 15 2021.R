
#loadin required packages 

library(dplyr)
library(tidyverse)
library(BiocManager)
library(edgeR)
library(Glimma)
library(tidyverse)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
library(emmeans)
library(qvalue)


###Count table import and preperation###

setwd("") #set working directory for the files


#read in count file
bluegill.full <- read_delim("Bluegill_subsetSuperTrans_counts_transcript_full.txt.gz",skip=1,delim="\t") 


view(table(duplicated(bluegill.full[,1]))) #Check for duplicate values in gene ID's

check2 <- bluegill.full %>% #another check for duplicate values
  group_by(Geneid) %>% 
  mutate(dupe = n()>1)

check3 <-str_detect(check2$dupe, "TRUE", negate = FALSE) %>% print() #goes with the above


#A final check for duplicate gene ID's

length (bluegill.full$Geneid) #109702 rows found
length(unique(bluegill.full$Geneid)) #109702 unique rows there no duplicate gene IDs' here
bluegill.full [duplicated(bluegill.full$Geneid), ] #0 duplicate values


#Extract label names from the text file

labels <- colnames(bluegill.full, do.NULL = TRUE, prefix = "col") 
write.csv(labels, "blugill labelsraw.csv") #write it as CSV to be able to extract the sample ID from the .bam file name


head(bluegill.full)
names(bluegill.full)

#As below, you can see that the .bam name has the sequencing ID but also the sample details which need to be isolated

# [1] "Geneid"                                                                                                            
# [2] "Chr"                                                                                                               
# [3] "Start"                                                                                                             
# [4] "End"                                                                                                               
# [5] "Strand"                                                                                                            
# [6] "Length"                                                                                                            
# [7] "mapped_trim_NS.1271.003.NEBNext_dual_i7_186---NEBNext_dual_i5_186.BC0-3G_R1.fastq.gzAligned.sortedByCoord.out.bam" 
# [8] "mapped_trim_NS.1271.003.NEBNext_dual_i7_187---NEBNext_dual_i5_187.BC0-4G_R1.fastq.gzAligned.sortedByCoord.out.bam" 
# [9] "mapped_trim_NS.1271.003.NEBNext_dual_i7_188---NEBNext_dual_i5_188.BC0-5G_R1.fastq.gzAligned.sortedByCoord.out.bam" 
# [10] "mapped_trim_NS.1271.003.NEBNext_dual_i7_189---NEBNext_dual_i5_189.BC0-8G_R1.fastq.gzAligned.sortedByCoord.out.bam" 
# [11] "mapped_trim_NS.1271.003.NEBNext_dual_i7_190---NEBNext_dual_i5_190.BC6-1G_R1.fastq.gzAligned.sortedByCoord.out.bam" 
# Geneid Chr   Start End   Strand Length `mapped_trim_NS~ `mapped_trim_NS~ `mapped_trim_NS~ `mapped_trim_NS~ `mapped_trim_NS~
#   <chr>  <chr> <chr> <chr> <chr>   <dbl>            <dbl>            <dbl>            <dbl>            <dbl>            <dbl>
#   1 Clust~ Clus~ 1;300 299;~ .;.       670               14               27              111               37              180
# 2 Clust~ Clus~ 1;27~ 277;~ .;.;.~   3302              444              648              972              474             1467
# 3 Clust~ Clus~ 1     655   .         655                1                0                0                0                0
# 4 Clust~ Clus~ 1;16~ 163;~ .;.;.~   2417               15               13               24                1               67
# 5 Clust~ Clus~ 1     1354  .        1354                1                3                6                3                7
# 6 Clust~ Clus~ 1;10~ 1033~ .;.;.~   2687               33              110               75               79              105





#Renaming the columns as the correct sample ID's 



bluegill.fullRN<-bluegill.full %>% rename('C0-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_186---NEBNext_dual_i5_186.BC0-3G_R1.fastq.gzAligned.sortedByCoord.out.bam', 
                                          'C0-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_187---NEBNext_dual_i5_187.BC0-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_188---NEBNext_dual_i5_188.BC0-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-8G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_189---NEBNext_dual_i5_189.BC0-8G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_190---NEBNext_dual_i5_190.BC6-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_191---NEBNext_dual_i5_191.BC6-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_192---NEBNext_dual_i5_192.BC6-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-7G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A11---NEBNext_dual_i5_A11.BC6-7G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-8G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.BC12-8G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.BT6-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A2---NEBNext_dual_i5_A2.BT12-3G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A3---NEBNext_dual_i5_A3.BT24-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A4---NEBNext_dual_i5_A4.BN6-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A5---NEBNext_dual_i5_A5.BN24-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A6---NEBNext_dual_i5_A6.BTN12-3G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-2L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A7---NEBNext_dual_i5_A7.BC0-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-2L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_A8---NEBNext_dual_i5_A8.BC6-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-8G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B11---NEBNext_dual_i5_B11.BC6-8G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B12---NEBNext_dual_i5_B12.BC24-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B1---NEBNext_dual_i5_B1.BT6-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B2---NEBNext_dual_i5_B2.BT12-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-7G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B3---NEBNext_dual_i5_B3.BT24-7G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B4---NEBNext_dual_i5_B4.BN6-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B5---NEBNext_dual_i5_B5.BTN6-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B6---NEBNext_dual_i5_B6.BTN12-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-3L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B7---NEBNext_dual_i5_B7.BC0-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-3L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_B8---NEBNext_dual_i5_B8.BC6-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C11---NEBNext_dual_i5_C11.BC12-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C12---NEBNext_dual_i5_C12.BC24-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C1---NEBNext_dual_i5_C1.BT6-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C2---NEBNext_dual_i5_C2.BT12-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-8G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C3---NEBNext_dual_i5_C3.BT24-8G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C4---NEBNext_dual_i5_C4.BN12-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C5---NEBNext_dual_i5_C5.BTN6-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C6---NEBNext_dual_i5_C6.BTN24-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-4L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_C7---NEBNext_dual_i5_C7.BC0-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D11---NEBNext_dual_i5_D11.BC12-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.BC24-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D1---NEBNext_dual_i5_D1.BT6-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D2---NEBNext_dual_i5_D2.BT12-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-9G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D3---NEBNext_dual_i5_D3.BT24-9G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D4---NEBNext_dual_i5_D4.BN12-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D5---NEBNext_dual_i5_D5.BTN6-3G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D6---NEBNext_dual_i5_D6.BTN24-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-5L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_D7---NEBNext_dual_i5_D7.BC0-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E11---NEBNext_dual_i5_E11.BC12-3G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E12---NEBNext_dual_i5_E12.BC24-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-7G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E1---NEBNext_dual_i5_E1.BT6-7G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-7G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E2---NEBNext_dual_i5_E2.BT12-7G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E3---NEBNext_dual_i5_E3.BN6-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E4---NEBNext_dual_i5_E4.BN12-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E5---NEBNext_dual_i5_E5.BTN6-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E6---NEBNext_dual_i5_E6.BTN24-3G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-6L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_E7---NEBNext_dual_i5_E7.BC0-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F11---NEBNext_dual_i5_F11.BC12-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F12---NEBNext_dual_i5_F12.BC24-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-9G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F1---NEBNext_dual_i5_F1.BT6-9G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-8G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F2---NEBNext_dual_i5_F2.BT12-8G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F3---NEBNext_dual_i5_F3.BN6-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F4---NEBNext_dual_i5_F4.BN12-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F5---NEBNext_dual_i5_F5.BTN6-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F6---NEBNext_dual_i5_F6.BTN24-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-7L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_F7---NEBNext_dual_i5_F7.BC0-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-5G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.BC12-5G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-7G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.BC24-7G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G1---NEBNext_dual_i5_G1.BT12-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-9G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G2---NEBNext_dual_i5_G2.BT12-9G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-3G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G3---NEBNext_dual_i5_G3.BN6-3G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G4---NEBNext_dual_i5_G4.BN24-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G5---NEBNext_dual_i5_G5.BTN12-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G6---NEBNext_dual_i5_G6.BTN24-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-8L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_G7---NEBNext_dual_i5_G7.BC0-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-1G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H10---NEBNext_dual_i5_H10.BC0-1G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-6G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H11---NEBNext_dual_i5_H11.BC12-6G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-8G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H12---NEBNext_dual_i5_H12.BC24-8G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H1---NEBNext_dual_i5_H1.BT12-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H2---NEBNext_dual_i5_H2.BT24-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H3---NEBNext_dual_i5_H3.BN6-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-4G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H4---NEBNext_dual_i5_H4.BN24-4G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-2G'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H5---NEBNext_dual_i5_H5.BTN12-2G_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C0-1L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H6---NEBNext_dual_i5_H6.BC0-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-1L'	=	'mapped_trim_NS.1271.003.NEBNext_dual_i7_H7---NEBNext_dual_i5_H7.BC6-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-8L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A10---NEBNext_dual_i5_A10.BC12-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-8L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A11---NEBNext_dual_i5_A11.BC24-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-9L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A12---NEBNext_dual_i5_A12.BT6-9L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-8L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.BT12-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-7L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A2---NEBNext_dual_i5_A2.BT24-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A3---NEBNext_dual_i5_A3.BN6-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A4---NEBNext_dual_i5_A4.BN24-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A5---NEBNext_dual_i5_A5.BTN6-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A6---NEBNext_dual_i5_A6.BTN24-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_A9---NEBNext_dual_i5_A9.BC12-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B10---NEBNext_dual_i5_B10.BC24-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B11---NEBNext_dual_i5_B11.BT6-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B12---NEBNext_dual_i5_B12.BT12-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-9L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B1---NEBNext_dual_i5_B1.BT12-9L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-8L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B2---NEBNext_dual_i5_B2.BT24-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B3---NEBNext_dual_i5_B3.BN12-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B4---NEBNext_dual_i5_B4.BN24-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B5---NEBNext_dual_i5_B5.BTN6-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B6---NEBNext_dual_i5_B6.BTN24-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_B9---NEBNext_dual_i5_B9.BC12-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C10---NEBNext_dual_i5_C10.BC24-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C11---NEBNext_dual_i5_C11.BT6-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C12---NEBNext_dual_i5_C12.BT12-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C1---NEBNext_dual_i5_C1.BT24-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-9L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C2---NEBNext_dual_i5_C2.BT24-9L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C3---NEBNext_dual_i5_C3.BN12-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C4---NEBNext_dual_i5_C4.BN24-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C5---NEBNext_dual_i5_C5.BTN12-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C6---NEBNext_dual_i5_C6.BTN24-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C8---NEBNext_dual_i5_C8.BC6-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_C9---NEBNext_dual_i5_C9.BC12-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D10---NEBNext_dual_i5_D10.BC24-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D11---NEBNext_dual_i5_D11.BT6-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D12---NEBNext_dual_i5_D12.BT12-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D1---NEBNext_dual_i5_D1.BT24-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D2---NEBNext_dual_i5_D2.BN6-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D3---NEBNext_dual_i5_D3.BN12-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D4---NEBNext_dual_i5_D4.BN24-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D5---NEBNext_dual_i5_D5.BTN12-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D6---NEBNext_dual_i5_D6.BTN24-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D8---NEBNext_dual_i5_D8.BC6-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_D9---NEBNext_dual_i5_D9.BC12-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E10---NEBNext_dual_i5_E10.BC24-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E11---NEBNext_dual_i5_E11.BT6-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E12---NEBNext_dual_i5_E12.BT12-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E1---NEBNext_dual_i5_E1.BT24-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E2---NEBNext_dual_i5_E2.BN6-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E3---NEBNext_dual_i5_E3.BN12-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E4---NEBNext_dual_i5_E4.BTN6-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E5---NEBNext_dual_i5_E5.BTN12-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E6---NEBNext_dual_i5_E6.BTN24-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E8---NEBNext_dual_i5_E8.BC6-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_E9---NEBNext_dual_i5_E9.BC12-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F10---NEBNext_dual_i5_F10.BC24-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F11---NEBNext_dual_i5_F11.BT6-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F12---NEBNext_dual_i5_F12.BT12-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F1---NEBNext_dual_i5_F1.BT24-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F2---NEBNext_dual_i5_F2.BN6-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F3---NEBNext_dual_i5_F3.BN12-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-2L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F4---NEBNext_dual_i5_F4.BTN6-2L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F5---NEBNext_dual_i5_F5.BTN12-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-7L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F8---NEBNext_dual_i5_F8.BC6-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-7L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_F9---NEBNext_dual_i5_F9.BC12-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G10---NEBNext_dual_i5_G10.BC24-6L_R1.fastq.gzAligned.sortedByCoord.out.bam', 
                                          'T6-7L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G11---NEBNext_dual_i5_G11.BT6-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.BT12-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G1---NEBNext_dual_i5_G1.BT24-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G2---NEBNext_dual_i5_G2.BN6-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N12-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G3---NEBNext_dual_i5_G3.BN12-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-3L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G4---NEBNext_dual_i5_G4.BTN6-3L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN12-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G5---NEBNext_dual_i5_G5.BTN12-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C6-8L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_G8---NEBNext_dual_i5_G8.BC6-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C24-7L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H10---NEBNext_dual_i5_H10.BC24-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T6-8L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H11---NEBNext_dual_i5_H11.BT6-8L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T12-7L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H12---NEBNext_dual_i5_H12.BT12-7L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'T24-6L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H1---NEBNext_dual_i5_H1.BT24-6L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N6-5L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H2---NEBNext_dual_i5_H2.BN6-5L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'N24-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H3---NEBNext_dual_i5_H3.BN24-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN6-4L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H4---NEBNext_dual_i5_H4.BTN6-4L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'TN24-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H5---NEBNext_dual_i5_H5.BTN24-1L_R1.fastq.gzAligned.sortedByCoord.out.bam',
                                          'C12-1L'	=	'mapped_trim_NS.1273.001.NEBNext_dual_i7_H8---NEBNext_dual_i5_H8.BC12-1L_R1.fastq.gzAligned.sortedByCoord.out.bam')
head(bluegill.fullRN) 

# 'C' = control, T = TFM, TN = TFM:niclosamide mixture, N = niclosamide, L = Liver, G = gill

#Checking that the renaming was applied                                                                                                                            

# Geneid Chr   Start End   Strand Length `C0-3G` `C0-4G` `C0-5G` `C0-8G` `C6-1G` `C6-4G` `C6-5G` `C6-7G` `C12-8G` `T6-1G`
# <chr>  <chr> <chr> <chr> <chr>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
#   1 Clust~ Clus~ 1;300 299;~ .;.       670      14      27     111      37     180      49      57     142      119      52
# 2 Clust~ Clus~ 1;27~ 277;~ .;.;.~   3302     444     648     972     474    1467     516     512    1482     1065     694
# 3 Clust~ Clus~ 1     655   .         655       1       0       0       0       0       0       0       1        0       0
# 4 Clust~ Clus~ 1;16~ 163;~ .;.;.~   2417      15      13      24       1      67      15      12      69       68      17
# 5 Clust~ Clus~ 1     1354  .        1354       1       3       6       3       7       1       5      10        6       3
# 6 Clust~ Clus~ 1;10~ 1033~ .;.;.~   2687      33     110      75      79     105      47      32     211       52      23









###########################data frame prep###################################


#First cleaning up the data frame a bit by removing 'Cluster' off of the Cluster ID's 
# renaming it ' geeID' and renaming  the length column 

erp <- select(bluegill.fullRN, -Strand, -Start, -End, -Chr) %>% 
  separate(Geneid, into=c("Cluster", "geeID"), sep = "-")%>% mutate(Lea= Length) %>% 
  unite("geeID", Cluster, geeID, sep = "_")

View(erp)
erp$geeID

table(duplicated(erp$geeID)) #Checking for duplicates
sum(is.na(erp$geeID)) #checking for NA's

#filtering out just the gill values (this is for gill TFM analyses)

Gilldata <- select(erp, contains("G"), contains("Lea")) %>% print()

dim(Gilldata) #74, 3 of these are geeID, lea and length

#THis project had some 0h controls that were not included in the analyses and have to be remove as below:

Gilldata2 <-select(Gilldata, -contains("C0")) %>% print()

dim(Gilldata2) #69

#Next removing the niclosamide and mixture fish from the data frame as 
#this one strictly dealt with TFM's effects in the gill

Gill_TFM <- select(Gilldata2, -contains("N"), -contains ("TN")) %>% print()

dim(Gill_TFM)
str(Gill_TFM )

#Reordering the columns for an easier order and removing the fish  'C24-5G' as there were issues with the sample

Gill_TFM  <- select(Gill_TFM , geeID, Lea, everything(),-contains("C24-5G")) 

dim(Gill_TFM) #40
str(Gill_TFM)

#All numbers match up with metadata and within the text file itself


#Below, getting the labels out to extract the treatment group and exposure duration information.
#This part was all done in excel 

gillTFMlabels <- colnames(Gill_TFM , do.NULL = TRUE, prefix = "col") %>% print()

write.csv(gillTFMlabels, "gillTFMlablesraw_Jan 4 2021.csv")


#Importing in the grouping information 
group_vector <- read.csv("groups_tfm_gills jan 4 2021.csv")

#uniting the group and time columns into a single column as a single factor (e.g. C.6, TFM.12, etc.)

grouping_info_TFM <-group_vector %>% 
  unite("treat_group_TFM",group:time, sep = ".", remove = TRUE, na.rm = FALSE) %>% 
  print()

#setting the grouping column as a factor
grouping_info_TFM$treat_group_TFM <- as.factor(grouping_info_TFM$treat_group_TFM) 


str(grouping_info_TFM) 
summary(grouping_info_TFM)






##########Using edgeR for differential expression ############

#Below, this will outline using the package 'edgeR' to input the count data into a DGElist object.
#It also includes the transformation, filtering, and quality control steps needed to prep the data for 
#differential expression analsyes in the edgeR package. 

dim(Gill_TFM)
str(Gill_TFM)

#making the count data into a DGElist object for edgeR
TFM_gill_edgeR <- DGEList(counts= Gill_TFM[,3:40], genes = Gill_TFM[, 1:2], group = grouping_info_TFM$treat_group_TFM) %>% 
  print()

colnames (TFM_gill_edgeR)
dim(TFM_gill_edgeR) #109702 38


###################step 1: Transformations from raw scale#####################


#Here, converting raw counts to counts per million (CPM) and log CPM within edgeR

cpm_values <- cpm(TFM_gill_edgeR) %>% print() #CPM
lcpm_values <- cpm(TFM_gill_edgeR, log=TRUE) %>% print() #log CPM


summary(cpm_values)
summary(lcpm_values)

#Getting the mean and median library sizes for the data 
L <- mean(TFM_gill_edgeR$samples$lib.size) * 1e-6
M <- median(TFM_gill_edgeR$samples$lib.size) * 1e-6
c(L, M) 

#[1] 32.77136 27.98592

#########Step 2: Filtering the data####################



#Setting up a design array for conducting filtering 

designTFM <- model.matrix(~0 + TFM_gill_edgeR$samples$group)
colnames(designTFM) <- c(levels(TFM_gill_edgeR$samples$group))
designTFM

#Next, we used the 'filterbyexp' function to filter the data based on the design matrix specified above and
#will remove lowly expressed genes

Keep.TFM <- filterByExpr(TFM_gill_edgeR, designTFM)

filt_TFM <- TFM_gill_edgeR[Keep.TFM,, keep.lib.sizes=FALSE] # the latter term as 'F' is needed to recompute library sizes

dim(filt_TFM)

dim(TFM_gill_edgeR) 



#########Step 3: Normalization of data#########################

#Post-filtering, the data is normalized using the  trimmed mean of M-values (TMM) 

normfact_TFM <- calcNormFactors(filt_TFM, method = "TMM") 

normfact_TFM$samples$norm.factors #looking at normalization factors just as a check


####Step 4: Data visualization####

#starting with a PCA to see how the data is grouping together 

#first, get CPM values for the plot
cpm_norm <- cpm(normfact_TFM) #getting CPM values
rownames(cpm_norm)<-normfact_TFM$genes$geeID
class(cpm_norm)
head(cpm_norm)

#Also extracting out the log CPM values for use in the heatmaps later on 

Log_CPM <-cpm(normfact_TFM, log =T)

normfact_TFM$genes$geeID

rownames(Log_CPM)<-normfact_TFM$genes$geeID
Log_CPM

Log_CPM.1 <- data.frame(Cluster = row.names(Log_CPM), Log_CPM)
head(Log_CPM.1)

write.csv(Log_CPM.1, "Log CPM_norm values_BG_TFM_gill.csv") #writing out CPM values for heatmaps later



#Back to PCA's
#Using CPM matrix to run a PCA with the package 'FactoMiner'

gill_PCA <-PCA(t(cpm_norm), graph = F) 

head(get_eig(gill_PCA)) #Getting the Eigenvalues

# eigenvalue variance.percent cumulative.variance.percent
# Dim.1   4356.314         9.436809                    9.436809
# Dim.2   3217.135         6.969077                   16.405886
# Dim.3   2921.712         6.329121                   22.735007
# Dim.4   2600.311         5.632890                   28.367897
# Dim.5   2106.839         4.563912                   32.931809
# Dim.6   2039.402         4.417827                   37.349636

#Visualizing the eigenvalues/variances

fviz_screeplot(gill_PCA, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = 10)

#Extracting the results for variables
var <- get_pca_var(gill_PCA) %>% print()

#Extracting the results for individuals
ind <- get_pca_ind(gill_PCA) %>% print()


#Actually making the PCA plot is below

#setting font to Times New Roman 
windowsFonts(Times=windowsFont("TT Times New Roman")) 

#renaming the 'normfact' object to avoid issues
New_dataF <-normfact_TFM 

#Setting up the group ID's for the plot and reording the groups 
New_dataF$samples$group <- factor(New_dataF $samples$group, levels = c("C.6", "C.12", "C.24", "T.6", "T.12", "T.24"))

#Plot for principle components 1 and 2

pca_PC1_PC2 <- fviz_pca_ind(gill_PCA,
                            geom.ind = "point", # show points only (but not "text"),
                            pointsize = 2.25,
                            legend.title = "", 
                            fill.ind = New_dataF$samples$group, # color by groups
                            col.ind = New_dataF$samples$group,
                            mean.point = F,  # removes group mean point
                            title = "",
                            addEllipses = F) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Times")) +
  scale_shape_manual(values = c(21, 22, 24, 21, 22, 24)) +
 theme(legend.position =  c(0.9, 0.25))+
   scale_fill_manual(values = c("blue1", "mediumpurple1",
                               "turquoise1", "gold1",
                               "orangered2", "darkorange3")) +
  scale_colour_manual(values = c("blue1", "mediumpurple1",
                               "turquoise1", "gold1",
                                 "orangered2", "darkorange3"))+
   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pca_PC1_PC2
ggsave("gill_TFM_PCA_Jan420201_simple.tiff", 
       plot = pca_PC1_PC2,
       width = 6, height = 6, units = "in", dpi = 300)




#MDS plotting

#Using glimma, can get more in depth by making MDS plots via glMDSPlot

glMDSPlot(normfact_TFM, groups=normfact_TFM$samples$group, html = "glMDS") 






######Step 5: Dispersion estimation#####

#Estimating dispersion in the data for differtianl expression analyses using the 'estimateDisp' function in edgeR

normfact_TFM <- estimateDisp(normfact_TFM, designTFM, robust = TRUE) %>% print()


#looking at some of the resulting dispersion parameters including 
#plotting the Biological Coefficient of Variation (BCV)

dim(designTFM) #38 6 
normfact_TFM$common.dispersion #0.1901501
summary(normfact_TFM$tagwise.dispersion)
normfact_TFM$prior.df

plotBCV(normfact_TFM)



#Next, taking the above and fitting it in a qusaliklihood glm

TFMfit <- glmQLFit(normfact_TFM, designTFM, robust=TRUE)


head(TFMfit$coefficients)
dim(TFMfit) #46163 6

#resulting coefficients from the glm

# C.12      C.24       C.6      T.12      T.24       T.6
# 1 -12.66830 -12.55320 -12.73225 -13.00937 -13.01553 -12.91125
# 2 -10.49313 -10.59787 -10.41646 -10.42089 -10.40679 -10.47755
# 4 -13.81832 -14.42284 -13.77129 -14.69689 -14.39411 -14.02912
# 6 -13.06317 -12.94567 -12.86754 -12.40302 -13.19899 -12.87931
# 8 -14.02287 -14.08619 -13.99255 -13.98079 -14.05011 -14.39910
# 9 -10.33556 -10.38589 -10.29130 -10.37249 -10.28859 -10.41736


#Can also visually inspect the result of the model as below:
plotQLDisp(TFMfit)

############Step 6: Detecting differential expression ######################################

#First, need to specify the contrasts of interest for the DE detection
#For this project, we're interested in control vs TFM at a given timepoint. 
#Making a matrix of our comparisons

cont.matrix_TFM <- makeContrasts(C.6vT.6 = T.6-C.6,     #control vs TFM @ 6h
                            C.12vT.12 = T.12 - C.12,    #control vs TFM @ 12h
                             C.24vT.24 = T.24 - C.24,   #Control vs TFM @ 24h
                            levels=designTFM) %>% print()


#Using an F test to make the comparisons for each of the contrasts and saving it as an object

CC.6vT.6 <- glmQLFTest(TFMfit, contrast = cont.matrix_TFM[ , 1])
CC.12vT.12 <- glmQLFTest(TFMfit, contrast = cont.matrix_TFM[ , 2])
CC.24vT.24 <- glmQLFTest(TFMfit, contrast = cont.matrix_TFM[ , 3])

#To get raw numbers of DE genes, coding it simply as 1, 0, -1 (up, null, down)

edgeR.coded_0.05_TFM<- new("TestResults", cbind(C.6vT.6 = decideTestsDGE(CC.6vT.6)[ , 1],
                                                C.12vT.12  = decideTestsDGE(CC.12vT.12)[ , 1],
                                                C.24vT.24 = decideTestsDGE(CC.24vT.24)[ , 1]))

# to keep the gene identifiers with the new object
rownames(edgeR.coded_0.05_TFM) <- rownames(CC.6vT.6$table) 

summary(edgeR.coded_0.05_TFM)

#Provides an out put of numbers of DE genes
# C.6vT.6 C.12vT.12 C.24vT.24
# -1     219      1061      2161
# 0    45568     43501     41793
# 1      376      1601      2209

#Can quickly visualize the results using the following plot for each comparison 

glMDPlot(CC.6vT.6, status=edgeR.coded_0.05_TFM[,1], counts= filt_TFM,
         groups=filt_TFM$samples$group, transform = TRUE, html = "CC.6vT.6")

glMDPlot(CC.12vT.12, status=edgeR.coded_0.05_TFM[,1], counts= filt_TFM,
         groups=filt_TFM$samples$group, transform = TRUE, html = "CC.12vT.12")

glMDPlot(CC.24vT.24, status=edgeR.coded_0.05_TFM[,1], counts= filt_TFM,
         groups=filt_TFM$samples$group, transform = TRUE, html = "CC.24vT.24")


#Getting detailed outputs for the contrasts
#Also includes calculating the fold change from the log fold change

# 1. Detailed results for the C6vT6 comparison and output to a file
CC.6vT.6.detailed <- topTags(CC.6vT.6, n = Inf, sort.by = "none")$table
CC.6vT.6.detailed$FC <- 2^abs(CC.6vT.6.detailed$logFC) * sign(CC.6vT.6.detailed$logFC)

CC.6vT.6.detailed

# 2. Detailed results for the C12vT12 comparison and output to a file
CC.12vT.12.detailed <- topTags(CC.12vT.12, n = Inf, sort.by = "none")$table
CC.12vT.12.detailed$FC <- 2^abs(CC.12vT.12.detailed$logFC) * sign(CC.12vT.12.detailed$logFC)

CC.12vT.12.detailed

# 3. Detailed results for the C24vT24 comparison and output to a file
CC.24vT.24.detailed <- topTags(CC.24vT.24, n = Inf, sort.by = "none")$table
CC.24vT.24.detailed$FC <- 2^abs(CC.24vT.24.detailed$logFC) * sign(CC.24vT.24.detailed$logFC)

CC.24vT.24.detailed


## Combining all of the detailed results along into one object that can be added to and then output and opened in Excel.

#Checking  to make sure they are in same order:
all.equal(rownames(CC.6vT.6.detailed), rownames(CC.12vT.12.detailed), rownames(CC.24vT.24.detailed))
# TRUE

all.equal(CC.6vT.6.detailed$geeID, CC.12vT.12.detailed$geeID, CC.24vT.24.detailed$geeID)


#Merge all of the results together and renaming the columns to specify the comparison:

head(CC.6vT.6.detailed[,c(1,4,3,8,7)])
head(CC.6vT.6.detailed)
dim(CC.6vT.6.detailed)
all.results_TFM <- bind_cols(gene_id = rownames(CC.6vT.6.detailed),
                             CC.6vT.6.detailed[,c(1,4,3,8,6, 7)], 
                             CC.12vT.12.detailed[,c(1,4,3,8,6, 7)], 
                             CC.24vT.24.detailed[,c(1,4,3,8,6, 7)]) %>% 
                             print()

  
#Consolidating all of the geeID's into single column and renaming the columns

all.results_TFM <- all.results_TFM %>% select(-gene_id, -geeID...8, -geeID...14) %>% #removing extra columns
   rename_at(1, funs(paste0("Cluster")))%>% #renaming geeID as Cluster for annotation joining
  rename_at(2:6, funs(paste0("C.6vT.6_", c("logCPM", "logFC", "FC", "PValue", "FDR")))) %>%  #renaming the columns to delineate by contrast
  rename_at(7:11, funs(paste0("C.12vT.12_", c("logCPM", "logFC", "FC", "PValue", "FDR")))) %>% 
  rename_at(12:16, funs(paste0("C.24vT.24_", c("logCPM", "logFC", "FC", "PValue", "FDR")))) %>% 
  print()

str(all.results_TFM)
names(all.results_TFM)
dim(all.results_TFM)  #46163 16, looks good here!

#Write out file so that it can be imported later
write_csv(all.results_TFM, "all.results_TFM_BG_gill.csv")


#####Step 7: Read in transcriptome and combine with DE results#######

#reading in a  copy of the annotated transcriptome
BGtranscriptome <- read.csv("bluegill_transcriptome.csv", header = T) 
dim(BGtranscriptome) #109702 41
  
 #keeping only clusters common to both results and transcriptome

BGtranscriptome <- BGtranscriptome %>% filter(Cluster %in% all.results_TFM$Cluster) %>%  
  print()

dim(BGtranscriptome) #46163 41, same as the all results file so its keep only genes that are showing up there
View(BGtranscriptome)
names(BGtranscriptome )


#Adding annotation info to the all.results object
#For ease, I have the all.results saved as an Excel that be loaded in quickly

all.results_TFM <- read_csv("C:/Users/m_law/Documents/work/U Manitoba/lampricide molecular project/Count table analyses/Bluegill analyses/TFM gills/Jan 4 20210 reanalysis/all.results_TFM_BG_gill.csv") %>% 
  print()

names(all.results_TFM)
#[1] "Cluster"          "C.6vT.6_logCPM"   "C.6vT.6_logFC"    "C.6vT.6_FC"       "C.6vT.6_PValue"   "C.6vT.6_FDR"     
#[7] "C.12vT.12_logCPM" "C.12vT.12_logFC"  "C.12vT.12_FC"     "C.12vT.12_PValue" "C.12vT.12_FDR"    "C.24vT.24_logCPM"
#[13] "C.24vT.24_logFC"  "C.24vT.24_FC"     "C.24vT.24_PValue" "C.24vT.24_FDR"      

dim(all.results_TFM)
# 46163        16

table(all.results_TFM$Cluster %in% BGtranscriptome$Cluster) # are all the ID terms present in both tables?
#  TRUE 
# 46163, all of the annotated genes are found in the results file

str(all.results_TFM$Cluster)
str(BGtranscriptome$Cluster)

#Have to make 'cluster' as a character the below won't work with different data types

BGtranscriptome$Cluster <-as.character(BGtranscriptome$Cluster) 

#Are the ID terms in the same order?
all.equal(all.results_TFM$Cluster, BGtranscriptome$Cluster) 
#FALSE

BGtranscriptome$Cluster

#reordering the dataframe to ascending cluster order 
BGtranscriptome1 <- BGtranscriptome %>% arrange(Cluster)   
head(BGtranscriptome1$Cluster)

all.results_TFM1 <- all.results_TFM %>%  arrange(Cluster)
head(all.results_TFM1$Cluster)

#are the ID terms in the same order?
all.equal(all.results_TFM1$Cluster, BGtranscriptome1$Cluster) #
#TRUE

#So all ID's from BGtranscriptome are in all.results_TFM, and they are in the same order, 
#so can just add the annotation info to the end of the all.results dataframe

#joining the two dataframes together
all.results_annot <- all.results_TFM1 %>% 
  left_join(BGtranscriptome1, by = "Cluster") %>% 
  print()
dim(all.results_annot)
# 46163 56


View(all.results_annot)

#write out the annotated results file

write_csv(all.results_annot, "all.results_annot.csv")




###########Step 8: Sorting and filtering the DE genes #####

#For ease reading in the annotated DE results file
all.results_annot <-read_csv("all.results_annot.csv")

#Just 6h fish

all.results_annot
names(all.results_annot)
View(all.results_annot)


#Grabbing columns for just 6h results and the corresponding gene/protein identifiers (i.e sprot files)
Gill_6h_DE<- all.results_annot %>% 
  select(1,2:6, 29:37, 17:28)
Gill_6h_DE
dim(Gill_6h_DE)
names(Gill_6h_DE)
#46163       27

#removing the rows with no annotation/gene.names
Gill_6h_DE.annot <- Gill_6h_DE %>% 
  filter(Gill_6h_DE$sp.BX_gene.name != "NA") #for now, only removing gene transcripts that don't have a gene name (ie. sp.bx)
dim(Gill_6h_DE.annot)                           
Gill_6h_DE.annot
#20817 27

#Filtering out  the genes that responded to treatment (i.e., significant genes)

Gill_6h_DE_SIG <- Gill_6h_DE.annot %>%
  filter(C.6vT.6_FDR <0.05)

dim(Gill_6h_DE_SIG)
#439 27

#Separating out the up regulated and downregulated genes
#Upregulated
Gill_6h_DE_UP <- Gill_6h_DE_SIG %>% 
  filter(C.6vT.6_FC > 0)

dim(Gill_6h_DE_UP)

#269 27

#Downregulated

Gill_6h_DE_DOWN <- Gill_6h_DE_SIG %>% 
  filter(C.6vT.6_FC < 0)

dim(Gill_6h_DE_DOWN)
#170 27




#Just 12 h fish


Gill_12h_DE<- all.results_annot %>% 
  select(1,7:11,29:37, 17:28)
names(Gill_12h_DE)
dim(Gill_12h_DE)
#46163       27

Gill_12h_DE.annot <- Gill_12h_DE %>% 
  filter(Gill_12h_DE$sp.BX_gene.name != "NA") 
dim(Gill_12h_DE.annot)                           
Gill_12h_DE.annot
#20817 27


Gill_12h_DE_SIG <- Gill_12h_DE.annot %>%
  filter(C.12vT.12_FDR <0.05)

dim(Gill_12h_DE_SIG)
#1808 27

#Upregulated
Gill_12h_DE_UP <- Gill_12h_DE_SIG %>% 
  filter(C.12vT.12_FC > 0)

dim(Gill_12h_DE_UP)

#1052 27

#Downregulated
Gill_12h_DE_DOWN <- Gill_12h_DE_SIG %>% 
  filter(C.12vT.12_FC < 0)

dim(Gill_12h_DE_DOWN)
#756 27




#Just 24 h fish


Gill_24h_DE<- all.results_annot %>% 
  select(1,12:16, 29:37, 17:28)
names(Gill_24h_DE)
dim(Gill_24h_DE)
#46163       27

Gill_24h_DE.annot <- Gill_24h_DE %>% 
  filter(Gill_24h_DE$sp.BX_gene.name != "NA") 
dim(Gill_24h_DE.annot)                           
Gill_24h_DE.annot
# 20817 27


Gill_24h_DE_SIG <- Gill_24h_DE.annot %>%
  filter(C.24vT.24_FDR <0.05)

dim(Gill_24h_DE_SIG)
#3158 27

#Upregulated genes 
Gill_24h_DE_UP <- Gill_24h_DE_SIG %>% 
  filter(C.24vT.24_FC > 0)

dim(Gill_24h_DE_UP)

#1546 27

#Downregulated genes

Gill_24h_DE_DOWN <- Gill_24h_DE_SIG %>% 
  filter(C.24vT.24_FC < 0)

dim(Gill_24h_DE_DOWN)
#1612 27



############EnrichR analyses############

#This series of analyses are used to determine if there is enrichment of gene oncology terms associated with 
#the DE genes and uses the package 'enrichR'

library(enrichR)


#First step is to make simplified objects containing just the false discovery rates (FDR), fold changes (FC), and 
#the gene names for each set of up/downregulated genes at each of the time points (i.e., 6, 12, 24 h)

#6h fish



#upregulated 


Gill_6h_DE_UP.Fin <-Gill_6h_DE_UP %>% select(C.6vT.6_FDR, C.6vT.6_FC, sp.BX_gene.name) %>% print() 
                       
              

#writing up as a .csv for getting gene lists
Gill_6h_DE_UP.Fin %>%
  write_tsv("TFM_BGgill_6h_UP.csv", col_names = F) 
dim(Gill_6h_DE_UP.Fin)
#269



#Downregulated


Gill_6h_DE_DOWN.Fin <-Gill_6h_DE_DOWN %>% select(C.6vT.6_FDR, C.6vT.6_FC, sp.BX_gene.name) 

Gill_6h_DE_DOWN.Fin %>% 
  write_tsv("TFM_BGgill_6h_DOWN.csv", col_names = F) 

dim(Gill_6h_DE_DOWN.Fin)
#170




# 12h fish

#Upregulated


Gill_12h_DE_UP.Fin <-Gill_12h_DE_UP %>% select(C.12vT.12_FDR, C.12vT.12_FC, sp.BX_gene.name) 


Gill_12h_DE_UP.Fin %>%
  write_tsv("TFM_BGgill_12h_UP.csv", col_names = F) 

dim(Gill_12h_DE_UP.Fin)
#1052



#Downregulated


Gill_12h_DE_DOWN.Fin <-Gill_12h_DE_DOWN %>% select(C.12vT.12_FDR, C.12vT.12_FC, sp.BX_gene.name) 

Gill_12h_DE_DOWN.Fin  %>% write_tsv("TFM_BGgill_12h_DOWN.csv", col_names = F) 

dim(Gill_12h_DE_DOWN.Fin)
#756



# 24h fish

#Upregulated


Gill_24h_DE_UP.Fin <-Gill_24h_DE_UP %>% select(C.24vT.24_FDR, C.24vT.24_FC, sp.BX_gene.name) 

Gill_24h_DE_UP.Fin  %>% 
  write_tsv("C:/Users/m_law/Documents/work/U Manitoba/lampricide molecular project/Count table analyses/Bluegill analyses/TFM gills/Jan 4 20210 reanalysis/REVIGO/24 h/TFM_BGgill_24h_UP.csv", col_names = F) 

dim(Gill_24h_DE_UP.Fin)
#1546



#Downregulated

Gill_24h_DE_DOWN.Fin <-Gill_24h_DE_DOWN %>% select(C.24vT.24_FDR, C.24vT.24_FC, sp.BX_gene.name) 


#writing up csv for getting gene lists
Gill_24h_DE_DOWN.Fin  %>% 
  write_tsv("C:/Users/m_law/Documents/work/U Manitoba/lampricide molecular project/Count table analyses/Bluegill analyses/TFM gills/Jan 4 20210 reanalysis/REVIGO/24 h/TFM_BGgill_24h_DOWN.csv", col_names = F) 

Gill_24h_DE_DOWN.Fin
dim(Gill_24h_DE_DOWN.Fin)
#1612 



####Running enrichR######


# For enrichR, need to first select the databases that will be used for the analyses
#We elected to use three databases as seen below


#Find available databases from EnrichR
dbs <- listEnrichrDbs() 

dbs <- c("GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "GO_Biological_Process_2018")

#Next, the gene lists for each contrast set were entered into the enrichR database 

# 6h fish 

# 6 h upregulated fish 
Gill_6h_DE_UP.Fin$sp.BX_gene.name <-as.character(Gill_6h_DE_UP.Fin$sp.BX_gene.name)

TFM_BG_gill_6h_UP_EnrichR<- enrichr(Gill_6h_DE_UP.Fin$sp.BX_gene.name, dbs) #actual Enrichr database comparison 

#Create tibbles for each GO category
TFM_BG_gill_6h_UP_BP <- as_tibble(TFM_BG_gill_6h_UP_EnrichR[["GO_Biological_Process_2018"]]) %>% 
  print()
TFM_BG_gill_6h_UP_MF <- as_tibble(TFM_BG_gill_6h_UP_EnrichR[["GO_Molecular_Function_2018"]]) %>% 
  print()
TFM_BG_gill_6h_UP_CC <- as_tibble(TFM_BG_gill_6h_UP_EnrichR[["GO_Cellular_Component_2018"]]) %>% 
  print()


#6 h DOWNregulated fish 

Gill_6h_DE_DOWN.Fin$sp.BX_gene.name <-as.character(Gill_6h_DE_DOWN.Fin$sp.BX_gene.name)

TFM_BG_gill_6h_DOWN_EnrichR<- enrichr(Gill_6h_DE_DOWN.Fin$sp.BX_gene.name, dbs) #actual Enrichr database comparison 

# Create tibbles for each GO category
TFM_BG_gill_6h_DOWN_BP <- as_tibble(TFM_BG_gill_6h_DOWN_EnrichR[["GO_Biological_Process_2018"]]) %>% 
  print()
TFM_BG_gill_6h_DOWN_MF <- as_tibble(TFM_BG_gill_6h_DOWN_EnrichR[["GO_Molecular_Function_2018"]]) %>% 
  print()
TFM_BG_gill_6h_DOWN_CC <- as_tibble(TFM_BG_gill_6h_DOWN_EnrichR[["GO_Cellular_Component_2018"]]) %>% 
  print()


# 12h fish 

# 12 h upregulated fish 
Gill_12h_DE_UP.Fin$sp.BX_gene.name <- as.character(Gill_12h_DE_UP.Fin$sp.BX_gene.name)

TFM_BG_gill_12h_UP_EnrichR<- enrichr(Gill_12h_DE_UP.Fin$sp.BX_gene.name, dbs) #actual Enrichr database comparison 

# Create tibbles for each GO category
TFM_BG_gill_12h_UP_BP <- as_tibble(TFM_BG_gill_12h_UP_EnrichR[["GO_Biological_Process_2018"]]) %>% 
  print()
TFM_BG_gill_12h_UP_MF <- as_tibble(TFM_BG_gill_12h_UP_EnrichR[["GO_Molecular_Function_2018"]]) %>% 
  print()
TFM_BG_gill_12h_UP_CC <- as_tibble(TFM_BG_gill_12h_UP_EnrichR[["GO_Cellular_Component_2018"]]) %>% 
  print()


#12 h DOWNregulated fish 
Gill_12h_DE_DOWN.Fin$sp.BX_gene.name <- as.character(Gill_12h_DE_DOWN.Fin$sp.BX_gene.name)

TFM_BG_gill_12h_DOWN_EnrichR<- enrichr(Gill_12h_DE_DOWN.Fin$sp.BX_gene.name, dbs) #actual Enrichr database comparison 

# Create tibbles for each GO category
TFM_BG_gill_12h_DOWN_BP <- as_tibble(TFM_BG_gill_12h_DOWN_EnrichR[["GO_Biological_Process_2018"]]) %>% 
  print()
TFM_BG_gill_12h_DOWN_MF <- as_tibble(TFM_BG_gill_12h_DOWN_EnrichR[["GO_Molecular_Function_2018"]]) %>% 
  print()
TFM_BG_gill_12h_DOWN_CC <- as_tibble(TFM_BG_gill_12h_DOWN_EnrichR[["GO_Cellular_Component_2018"]]) %>% 
  print()


# 24h fish 

## 24 h upregulated fish 
Gill_24h_DE_UP.Fin$sp.BX_gene.name <-as.character(Gill_24h_DE_UP.Fin$sp.BX_gene.name)

TFM_BG_gill_24h_UP_EnrichR<- enrichr(Gill_24h_DE_UP.Fin$sp.BX_gene.name, dbs) #actual Enrichr database comparison 

# Create tibbles for each GO category
TFM_BG_gill_24h_UP_BP <- as_tibble(TFM_BG_gill_24h_UP_EnrichR[["GO_Biological_Process_2018"]]) %>% 
  print()
TFM_BG_gill_24h_UP_MF <- as_tibble(TFM_BG_gill_24h_UP_EnrichR[["GO_Molecular_Function_2018"]]) %>% 
  print()
TFM_BG_gill_24h_UP_CC <- as_tibble(TFM_BG_gill_24h_UP_EnrichR[["GO_Cellular_Component_2018"]]) %>% 
  print()


## 24 h DOWNregulated fish 
Gill_24h_DE_DOWN.Fin$sp.BX_gene.name <- as.character(Gill_24h_DE_DOWN.Fin$sp.BX_gene.name)
TFM_BG_gill_24h_DOWN_EnrichR<- enrichr(Gill_24h_DE_DOWN.Fin$sp.BX_gene.name, dbs) #actual Enrichr database comparison 

# Create tibbles for each GO category
TFM_BG_gill_24h_DOWN_BP <- as_tibble(TFM_BG_gill_24h_DOWN_EnrichR[["GO_Biological_Process_2018"]]) %>% 
  print()
TFM_BG_gill_24h_DOWN_MF <- as_tibble(TFM_BG_gill_24h_DOWN_EnrichR[["GO_Molecular_Function_2018"]]) %>% 
  print()
TFM_BG_gill_24h_DOWN_CC <- as_tibble(TFM_BG_gill_24h_DOWN_EnrichR[["GO_Cellular_Component_2018"]]) %>% 
  print()



########### Prepare files for Revigo ######

#Once the enrichR database comparison was made, the files needed to be prepped for entry into Revigo, which we are 
#using as the primary means of assaying GO term enrichment 

#In the below, this is simply organizing the dataframes for each timepoint/erichR data series to be read by revigo and
#to filter out significant GO terms with >= 4 genes. This was repeated for 6h, 12h, and 24h fish


#6 h TFM upregulated genes 

# Split Term into GO.Desc and GO.Term, add column for Gene.no
str (TFM_BG_gill_6h_UP_BP) #1804
str(TFM_BG_gill_6h_UP_MF) #335
str(TFM_BG_gill_6h_UP_CC) #109

## check that each term can be separated by "(GO:"

 
strsplit(TFM_BG_gill_6h_UP_BP$Term, " \\(GO:") %>% sapply(length) %>% table  #1804

strsplit(TFM_BG_gill_6h_UP_MF$Term, " \\(GO:") %>% sapply(length) %>% table #335

strsplit(TFM_BG_gill_6h_UP_CC$Term, " \\(GO:") %>% sapply(length) %>% table #109

TFM_BG_gill_6h_UP_BP.2 <- TFM_BG_gill_6h_UP_BP %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_6h_UP_BP$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_6h_UP_MF.2 <- TFM_BG_gill_6h_UP_MF %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_6h_UP_MF$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_6h_UP_CC.2 <- TFM_BG_gill_6h_UP_CC %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_6h_UP_CC$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()


# Create a table with just the necessary columns and significant terms (p < 0.5) and with >= 4 genes for REVIGO analysis
TFM_BG_gill_6h_UP_BP.3 <- TFM_BG_gill_6h_UP_BP.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_6h_UP_BP.3)
#43   5


TFM_BG_gill_6h_UP_MF.3 <- TFM_BG_gill_6h_UP_MF.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_6h_UP_MF.3)
#24   5

TFM_BG_gill_6h_UP_CC.3 <- TFM_BG_gill_6h_UP_CC.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_6h_UP_CC.3)
#0   5

# bind tables together and write out
bind_rows(TFM_BG_gill_6h_UP_BP.3, TFM_BG_gill_6h_UP_MF.3, TFM_BG_gill_6h_UP_CC.3) %>% 
  select(GO.Term, Adjusted.P.value) %>% 
  write_csv("REVIGO_BG_Gill_TFM_6h_UP.csv")


#6 h TFM DOWNregulated genes 

# Split Term into GO.Desc and GO.Term, add column for Gene.no
str (TFM_BG_gill_6h_DOWN_BP) #1291
str(TFM_BG_gill_6h_DOWN_MF) #215
str(TFM_BG_gill_6h_DOWN_CC) #83

## check that each term can be separated by "(GO:"


strsplit(TFM_BG_gill_6h_DOWN_BP$Term, " \\(GO:") %>% sapply(length) %>% table  #1291

strsplit(TFM_BG_gill_6h_DOWN_MF$Term, " \\(GO:") %>% sapply(length) %>% table #215

strsplit(TFM_BG_gill_6h_DOWN_CC$Term, " \\(GO:") %>% sapply(length) %>% table #83

TFM_BG_gill_6h_DOWN_BP.2 <- TFM_BG_gill_6h_DOWN_BP %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_6h_DOWN_BP$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_6h_DOWN_MF.2 <- TFM_BG_gill_6h_DOWN_MF %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_6h_DOWN_MF$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_6h_DOWN_CC.2 <- TFM_BG_gill_6h_DOWN_CC %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_6h_DOWN_CC$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()


# Create a table with just the necessary columns and significant terms (p < 0.5) and with >= 4 genes for REVIGO analysis
TFM_BG_gill_6h_DOWN_BP.3 <- TFM_BG_gill_6h_DOWN_BP.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_6h_DOWN_BP.3)
#44   5


TFM_BG_gill_6h_DOWN_MF.3 <- TFM_BG_gill_6h_DOWN_MF.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_6h_DOWN_MF.3)
#5   5

TFM_BG_gill_6h_DOWN_CC.3 <- TFM_BG_gill_6h_DOWN_CC.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_6h_DOWN_CC.3)
#0   5

# bind tables together and write out
bind_rows(TFM_BG_gill_6h_DOWN_BP.3, TFM_BG_gill_6h_DOWN_MF.3, TFM_BG_gill_6h_DOWN_CC.3) %>% 
  select(GO.Term, Adjusted.P.value) %>% 
  write_csv("REVIGO_BG_Gill_TFM_6h_DOWN.csv")





#12 h TFM upregulated genes 

# Split Term into GO.Desc and GO.Term, add column for Gene.no
str (TFM_BG_gill_12h_UP_BP) #3304
str(TFM_BG_gill_12h_UP_MF) #652
str(TFM_BG_gill_12h_UP_CC) #223

## check that each term can be separated by "(GO:"


strsplit(TFM_BG_gill_12h_UP_BP$Term, " \\(GO:") %>% sapply(length) %>% table  #3304

strsplit(TFM_BG_gill_12h_UP_MF$Term, " \\(GO:") %>% sapply(length) %>% table #652

strsplit(TFM_BG_gill_12h_UP_CC$Term, " \\(GO:") %>% sapply(length) %>% table #223

TFM_BG_gill_12h_UP_BP.2 <- TFM_BG_gill_12h_UP_BP %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_12h_UP_BP$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_12h_UP_MF.2 <- TFM_BG_gill_12h_UP_MF %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_12h_UP_MF$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_12h_UP_CC.2 <- TFM_BG_gill_12h_UP_CC %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_12h_UP_CC$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()


# Create a table with just the necessary columns and significant terms (p < 0.5) and with >= 4 genes for REVIGO analysis
TFM_BG_gill_12h_UP_BP.3 <- TFM_BG_gill_12h_UP_BP.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_12h_UP_BP.3)
#159  5


TFM_BG_gill_12h_UP_MF.3 <- TFM_BG_gill_12h_UP_MF.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_12h_UP_MF.3)
#30   5

TFM_BG_gill_12h_UP_CC.3 <- TFM_BG_gill_12h_UP_CC.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_12h_UP_CC.3)
#8   5

# bind tables together and write out
  bind_rows(TFM_BG_gill_12h_UP_BP.3, TFM_BG_gill_12h_UP_MF.3, TFM_BG_gill_12h_UP_CC.3) %>% 
  select(GO.Term, Adjusted.P.value) %>% 
  write_csv("REVIGO_BG_Gill_TFM_12h_UP.csv")


# 12h TFM DOWNregulated genes 

# Split Term into GO.Desc and GO.Term, add column for Gene.no
str (TFM_BG_gill_12h_DOWN_BP) #2877
str(TFM_BG_gill_12h_DOWN_MF) #610
str(TFM_BG_gill_12h_DOWN_CC) #232

## check that each term can be separated by "(GO:"


strsplit(TFM_BG_gill_12h_DOWN_BP$Term, " \\(GO:") %>% sapply(length) %>% table  #2877

strsplit(TFM_BG_gill_12h_DOWN_MF$Term, " \\(GO:") %>% sapply(length) %>% table #610

strsplit(TFM_BG_gill_12h_DOWN_CC$Term, " \\(GO:") %>% sapply(length) %>% table #232

TFM_BG_gill_12h_DOWN_BP.2 <- TFM_BG_gill_12h_DOWN_BP %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_12h_DOWN_BP$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_12h_DOWN_MF.2 <- TFM_BG_gill_12h_DOWN_MF %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_12h_DOWN_MF$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_12h_DOWN_CC.2 <- TFM_BG_gill_12h_DOWN_CC %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_12h_DOWN_CC$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()


# Create a table with just the necessary columns and significant terms (p < 0.5) and with >= 4 genes for REVIGO analysis
TFM_BG_gill_12h_DOWN_BP.3 <- TFM_BG_gill_12h_DOWN_BP.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_12h_DOWN_BP.3)
#40   5


TFM_BG_gill_12h_DOWN_MF.3 <- TFM_BG_gill_12h_DOWN_MF.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_12h_DOWN_MF.3)
#0   5

TFM_BG_gill_12h_DOWN_CC.3 <- TFM_BG_gill_12h_DOWN_CC.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_12h_DOWN_CC.3)
#0   5

# bind tables together and write out
bind_rows(TFM_BG_gill_12h_DOWN_BP.3, TFM_BG_gill_12h_DOWN_MF.3, TFM_BG_gill_12h_DOWN_CC.3) %>% 
  select(GO.Term, Adjusted.P.value) %>% 
  write_csv("REVIGO_BG_Gill_TFM_12h_DOWN.csv")






# 24 h TFM upregulated genes 

# Split Term into GO.Desc and GO.Term, add column for Gene.no
str (TFM_BG_gill_24h_UP_BP) #3732
str(TFM_BG_gill_24h_UP_MF) #797
str(TFM_BG_gill_24h_UP_CC) #285

## check that each term can be separated by "(GO:"


strsplit(TFM_BG_gill_24h_UP_BP$Term, " \\(GO:") %>% sapply(length) %>% table  #3732

strsplit(TFM_BG_gill_24h_UP_MF$Term, " \\(GO:") %>% sapply(length) %>% table #797

strsplit(TFM_BG_gill_24h_UP_CC$Term, " \\(GO:") %>% sapply(length) %>% table #285

TFM_BG_gill_24h_UP_BP.2 <- TFM_BG_gill_24h_UP_BP %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_24h_UP_BP$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_24h_UP_MF.2 <- TFM_BG_gill_24h_UP_MF %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_24h_UP_MF$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_24h_UP_CC.2 <- TFM_BG_gill_24h_UP_CC %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_24h_UP_CC$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()


# Create a table with just the necessary columns and significant terms (p < 0.5) and with >= 4 genes for REVIGO analysis
TFM_BG_gill_24h_UP_BP.3 <- TFM_BG_gill_24h_UP_BP.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_24h_UP_BP.3)
#82   5


TFM_BG_gill_24h_UP_MF.3 <- TFM_BG_gill_24h_UP_MF.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_24h_UP_MF.3)
#11   5

TFM_BG_gill_24h_UP_CC.3 <- TFM_BG_gill_24h_UP_CC.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_24h_UP_CC.3)
#9  5

# bind tables together and write out
bind_rows(TFM_BG_gill_24h_UP_BP.3, TFM_BG_gill_24h_UP_MF.3, TFM_BG_gill_24h_UP_CC.3) %>% 
  select(GO.Term, Adjusted.P.value) %>% 
  write_csv("REVIGO_BG_Gill_TFM_24h_UP.csv")


### 24 h TFM DOWNregulated genes ###

# Split Term into GO.Desc and GO.Term, add column for Gene.no
str (TFM_BG_gill_24h_DOWN_BP) #3927
str(TFM_BG_gill_24h_DOWN_MF) #813
str(TFM_BG_gill_24h_DOWN_CC) #338

## check that each term can be separated by "(GO:"


strsplit(TFM_BG_gill_24h_DOWN_BP$Term, " \\(GO:") %>% sapply(length) %>% table  #3927

strsplit(TFM_BG_gill_24h_DOWN_MF$Term, " \\(GO:") %>% sapply(length) %>% table #813

strsplit(TFM_BG_gill_24h_DOWN_CC$Term, " \\(GO:") %>% sapply(length) %>% table #338

TFM_BG_gill_24h_DOWN_BP.2 <- TFM_BG_gill_24h_DOWN_BP %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_24h_DOWN_BP$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_24h_DOWN_MF.2 <- TFM_BG_gill_24h_DOWN_MF %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_24h_DOWN_MF$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()

TFM_BG_gill_24h_DOWN_CC.2 <- TFM_BG_gill_24h_DOWN_CC %>%
  separate(Term, c("GO.Desc", "GO.Term"), " \\(GO:") %>%  # separate GO description and GO term
  separate(GO.Term, "GO.Term", "\\)", extra = "drop") %>% # remove extra ) at end of GO.Term
  mutate(GO.Term = paste0("GO:", GO.Term)) %>%   # add "GO:" back in front of GO.term
  mutate(Gene.no = strsplit(TFM_BG_gill_24h_DOWN_CC$Genes, ";") %>% sapply(length)) %>% # create a column with number of genes for each term
  print()


# Create a table with just the necessary columns and significant terms (p < 0.5) and with >= 4 genes for REVIGO analysis
TFM_BG_gill_24h_DOWN_BP.3 <- TFM_BG_gill_24h_DOWN_BP.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_24h_DOWN_BP.3)
#325   5


TFM_BG_gill_24h_DOWN_MF.3 <- TFM_BG_gill_24h_DOWN_MF.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_24h_DOWN_MF.3)
#52   5

TFM_BG_gill_24h_DOWN_CC.3 <- TFM_BG_gill_24h_DOWN_CC.2 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(Gene.no >= 4) %>% 
  dplyr::select(GO.Term, GO.Desc, Adjusted.P.value, Gene.no, Genes) %>% 
  print()
dim(TFM_BG_gill_24h_DOWN_CC.3)
#70   5

# bind tables together and write out
bind_rows(TFM_BG_gill_24h_DOWN_BP.3, TFM_BG_gill_24h_DOWN_MF.3, TFM_BG_gill_24h_DOWN_CC.3) %>% 
  select(GO.Term, Adjusted.P.value) %>% 
  write_csv("REVIGO_BG_Gill_TFM_24h_DOWN.csv")


#These files were then input into Revigo (http://revigo.irb.hr/) where the final part of the enrichment analysis occurred
#We used a similarity value of 0.5 (small)
#Below details how the plots were made for the resulting revigo summaries

#### Revigo plots ####

library(devtools)
library(treemapify)
library(reshape2)
library(tidyverse)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(edgeR)
library(Glimma)
library(tidyverse)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
library(emmeans)
library(qvalue)
library(ggplot2)
library(RColorBrewer)



setwd("")

#This outlines how the resulting Revigo data was made into the GO term count plots as in the manuscript. 
#For simplicity, we opted to restrict plot making to just the results for the biological processes (BP) and 
#molecular functions (MF). Cellular component (CC) was not used in the plots. 

#Figure data preperation

#This just involves curating the raw Revigo file into something that can be used in making the plots

#6 h

#UP


Revigo_BP_UP_6_05 <- read_csv("BG_gill_TFM_6 h_UP_BP_05_Jan 8.CSV", skip =4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Biological Processes") %>% 
  print()


Revigo_MF_UP_6_05 <- read_csv("BG_gill_TFM_6 h_UP_MF_05_Jan 8.CSV", skip =4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Molecular Functions") %>% 
  print()

#DOWN



Revigo_BP_DOWN_6_05 <- read_csv("BG_gill_TFM_6 h_DOWN_BP_05_Jan 8.CSV", skip =4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Biological Process") %>% 
  print()

Revigo_MF_DOWN_6_05 <- read_csv("BG_gill_TFM_6 h_DOWN_MF_05_Jan 8.CSV", skip =4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Molecular Functions") %>% 
  print()





# 12 h 


#UP


Revigo_BP_UP_12_05 <- read_csv("BG_gill_TFM_12 h_UP_BP_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Biological Processes") %>% 
  print()

Revigo_MF_UP_12_05 <- read_csv("BG_gill_TFM_12 h_UP_MF_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Molecular Functions") %>% 
  print()



#DOWN


Revigo_BP_DOWN_12_05 <- read_csv("BG_gill_TFM_12 h_DOWN_BP_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Biological Processes") %>% 
  print()

Revigo_MF_DOWN_12_05 <- read_csv("BG_gill_TFM_12 h_DOWN_MF_05_jan 8.csv", skip = 4) %>%
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>%
  mutate(frequency = as.numeric(frequency)) %>%
  mutate(abslog10pvalue = abs(log10pvalue)) %>%
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>%
  mutate(GO.cat = "Molecular Functions") %>%
  print()




#24 h 

#UP


Revigo_BP_UP_24_05 <- read_csv("BG_gill_TFM_24 h_UP_BP_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Biological Processes") %>% 
  print()

Revigo_MF_UP_24_05 <- read_csv("BG_gill_TFM_24 h_UP_MF_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Molecular Functions") %>% 
  print()



#DOWN



Revigo_BP_DOWN_24_05 <- read_csv("BG_gill_TFM_24 h_DOWN_BP_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Biological Processes") %>% 
  print()

Revigo_MF_DOWN_24_05 <- read_csv("BG_gill_TFM_24 h_DOWN_MF_05_jan 8.csv", skip = 4) %>% 
  separate(frequencyInDb, into = "frequency", sep = "%", extra = "drop") %>% 
  mutate(frequency = as.numeric(frequency)) %>% 
  mutate(abslog10pvalue = abs(log10pvalue)) %>% 
  select(term_ID:frequency, abslog10pvalue, uniqueness:representative) %>% 
  mutate(GO.cat = "Molecular Functions") %>% 
  print()




#####Creating summary tables for Revigo plots######

##creating summary tables of all of the above datframes for ease of reimport at a later time

#change col names so can join with the enrichR files

#6 h


Revigo_BP_UP_6_05.1 <- Revigo_BP_UP_6_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()
Revigo_MF_UP_6_05.1 <- Revigo_MF_UP_6_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()


Revigo_BP_DOWN_6_05.1 <- Revigo_BP_DOWN_6_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()
Revigo_MF_DOWN_6_05.1 <- Revigo_MF_DOWN_6_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()



#12 h


Revigo_BP_UP_12_05.1 <- Revigo_BP_UP_12_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()
Revigo_MF_UP_12_05.1 <- Revigo_MF_UP_12_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()



Revigo_BP_DOWN_12_05.1 <- Revigo_BP_DOWN_12_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()
Revigo_MF_DOWN_12_05.1 <- Revigo_MF_DOWN_12_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()


#24


Revigo_BP_UP_24_05.1 <- Revigo_BP_UP_24_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()
Revigo_MF_UP_24_05.1 <- Revigo_MF_UP_24_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()


Revigo_BP_DOWN_24_05.1 <- Revigo_BP_DOWN_24_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()
Revigo_MF_DOWN_24_05.1 <- Revigo_MF_DOWN_24_05 %>% select(GO.Term = term_ID, GO.Desc = description, GO.Rep = representative) %>% print()




#Joining the above objects with their corresponding enrichR files


#6 h


#Upregulated

Revigo_BP_UP_6_05_summary <- TFM_BG_gill_6h_UP_BP.3 %>%
  left_join(Revigo_BP_UP_6_05.1)  %>%
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>%
  mutate(GO.cat = "BP") %>%  #add row for GO category
  print()

Revigo_MF_UP_6_05_summary <- TFM_BG_gill_6h_UP_MF.3 %>%
  left_join(Revigo_MF_UP_6_05.1)  %>%
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>%
  mutate(GO.cat = "MF") %>%  #add row for GO category
  print()




#Downregulated
Revigo_BP_DOWN_6_05_summary <- TFM_BG_gill_6h_DOWN_BP.3 %>%
  left_join(Revigo_BP_DOWN_6_05.1)  %>%
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>%
  mutate(GO.cat = "BP") %>%  #add row for GO category
  print()

Revigo_MF_DOWN_6_05_summary <- TFM_BG_gill_6h_DOWN_MF.3 %>%
  left_join(Revigo_MF_DOWN_6_05.1)  %>%
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>%
  mutate(GO.cat = "MF") %>%  #add row for GO category
  print()





#12 h


#upregulated
Revigo_BP_UP_12_05_summary <- TFM_BG_gill_12h_UP_BP.3 %>% 
  left_join(Revigo_BP_UP_12_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "BP") %>% # add row for GO category
  print()

Revigo_MF_UP_12_05_summary <- TFM_BG_gill_12h_UP_MF.3 %>% 
  left_join(Revigo_MF_UP_12_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "MF") %>% # add row for GO category
  print() 



#Downregulated
Revigo_BP_DOWN_12_05_summary <- TFM_BG_gill_12h_DOWN_BP.3 %>% 
  left_join(Revigo_BP_DOWN_12_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "BP") %>% # add row for GO category
  print()

Revigo_MF_DOWN_12_05_summary <- TFM_BG_gill_12h_DOWN_MF.3 %>%
  left_join(Revigo_MF_DOWN_12_05.1)  %>%
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>%
  mutate(GO.cat = "MF") %>% # add row for GO category
  print()




#24 h

#Upregulated
Revigo_BP_UP_24_05_summary <- TFM_BG_gill_24h_UP_BP.3 %>% 
  left_join(Revigo_BP_UP_24_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "BP") %>% # add row for GO category
  print()

Revigo_MF_UP_24_05_summary <- TFM_BG_gill_24h_UP_MF.3 %>% 
  left_join(Revigo_MF_UP_24_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "MF") %>% # add row for GO category
  print() 


#Downregulated 
Revigo_BP_DOWN_24_05_summary <- TFM_BG_gill_24h_DOWN_BP.3 %>% 
  left_join(Revigo_BP_DOWN_24_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "BP") %>% # add row for GO category
  print()

Revigo_MF_DOWN_24_05_summary <- TFM_BG_gill_24h_DOWN_MF.3 %>% 
  left_join(Revigo_MF_DOWN_24_05.1)  %>% 
  select(GO.Rep, Genes, Gene.no, GO.Term:Adjusted.P.value) %>% 
  mutate(GO.cat = "MF") %>% # add row for GO category
  print() 




#Combine all of these dataframes and write out files as a .csv for ease




Revigo_BP_UP_6_05_summary %>% 
  bind_rows(Revigo_MF_UP_6_05_summary) %>% 
  write.csv("Gill_TFM_6_UP_GO_05_ALL_summary.csv")



Revigo_BP_DOWN_6_05_summary %>% 
  bind_rows(Revigo_MF_DOWN_6_05_summary) %>% 
  write.csv("Gill_TFM_6_DOWN_GO_05_ALL_summary.csv")



#12 h


Revigo_BP_UP_12_05_summary %>% 
  bind_rows(Revigo_MF_UP_12_05_summary) %>% 
  write.csv("Gill_TFM_12_UP_GO_05_ALL_summary.csv")


Revigo_BP_DOWN_12_05_summary %>% 
  bind_rows(Revigo_MF_DOWN_12_05_summary) %>% 
  write.csv("Gill_TFM_12_DOWN_GO_05_ALL_summary.csv")


#24 h


Revigo_BP_UP_24_05_summary %>% 
  bind_rows(Revigo_MF_UP_24_05_summary) %>% 
  write.csv("Gill_TFM_24_UP_GO_05_ALL_summary.csv")


Revigo_BP_DOWN_24_05_summary %>% 
  bind_rows(Revigo_MF_DOWN_24_05_summary) %>% 
  write.csv("Gill_TFM_24_DOWN_GO_05_ALL_summary.csv")




####Making the GO plots#####

#Read in the summary files made above

# 6 h



Gill_TFM_6_UP_GO_05_ALL_summary <- read.csv("Gill_TFM_6_UP_GO_05_ALL_summary.csv")

Gill_TFM_6_DOWN_GO_05_ALL_summary <-  read.csv("Gill_TFM_6_DOWN_GO_05_ALL_summary.csv")



#12 h

Gill_TFM_12_UP_GO_05_ALL_summary<-  read.csv("Gill_TFM_12_UP_GO_05_ALL_summary.csv")

Gill_TFM_12_DOWN_GO_05_ALL_summary<- read.csv("Gill_TFM_12_DOWN_GO_05_ALL_summary.csv")


#24 h


Gill_TFM_24_UP_GO_05_ALL_summary <-  read.csv("Gill_TFM_24_UP_GO_05_ALL_summary.csv")

Gill_TFM_24_DOWN_GO_05_ALL_summary <- read.csv("Gill_TFM_24_DOWN_GO_05_ALL_summary.csv")



#These functions below look to prep the above dataframes by removing GO terms that were removed by Revigo
#while re-leveling the major terms based on the representative GO term (GO.Rep). Also, adding in identifying 
#columns for fill colours and rearranging the tibble ina more favourable way.
#Also the series of mutate functions following the GO.Cat lines (i.e. past mutate(Colour = case_when(
#GO.cat %in% "BP" ~ "#99d8c9", GO.cat %in% "MF" ~ "#3690c0") are needed for the specific style of the journal 
#which wanted capital letters only on the first word (i.e. the 'str_to_sentence'). Then, in order to correct for
#words like 'DNA', 'GTP', 'ATP', etc, these were manually put in through the use of the 'str_replace' function which is 
#why the second half of preperation is so extensive




# 6 h 

#Upregualtion
Gill_TFM_6_UP_GO_05_ALL_summary.1 <- Gill_TFM_6_UP_GO_05_ALL_summary %>% 
  filter(!is.na(GO.Rep)) %>%  # remove terms eliminated by REVIGO
  mutate(GO.cat = factor(GO.cat, levels = c("BP", "MF")),
         GO.Rep = factor(GO.Rep, levels = c('cell junction assembly',
                                            'heart development',
                                            'regulation of transcription from RNA polymerase II promoter',
                                            'chemokine receptor activity',
                                            'phosphatidylinositol 3-kinase activity',
                                            'RNA polymerase II transcription factor binding',
                                            'transcription factor activity, RNA polymerase II distal enhancer sequence-specific binding',
                                            'transcription regulatory region DNA binding'
         )),
         GO.Desc = factor(GO.Desc)) %>% 
  arrange(GO.cat, GO.Rep, desc(Gene.no), GO.Desc) %>%
  select( GO.cat, GO.Rep, GO.Desc, GO.Term, Gene.no, Genes) %>% 
  mutate(Colour = case_when(
    GO.cat %in% "BP" ~ "#99d8c9",
    GO.cat %in% "MF" ~ "#3690c0"
  )) %>% mutate(GO.Rep = str_to_sentence(GO.Rep)) %>% 
  mutate(GO.Desc = str_to_sentence(GO.Desc)) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "dna", "DNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "ii", "II")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "ii", "II")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mapk", "MAPK")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mapk", "MAPK")) %>%
  print(width = Inf) 


#Downregualtion
Gill_TFM_6_DOWN_GO_05_ALL_summary.1 <- Gill_TFM_6_DOWN_GO_05_ALL_summary %>% 
  filter(!is.na(GO.Rep)) %>%  # remove terms eliminated by REVIGO
  mutate(GO.cat = factor(GO.cat, levels = c("BP", "MF")),
         GO.Rep = factor(GO.Rep, levels = c('actin cytoskeleton reorganization',
                                            'negative regulation of multicellular organismal process',
                                            'proteoglycan biosynthesis',
                                            'G-protein coupled receptor binding',
                                            'transmembrane receptor protein tyrosine kinase activity'
         )),
         GO.Desc = factor(GO.Desc)) %>% 
  arrange(GO.cat, GO.Rep, desc(Gene.no), GO.Desc) %>%
  select( GO.cat, GO.Rep, GO.Desc, GO.Term, Gene.no, Genes) %>% 
  mutate(Colour = case_when(
    GO.cat %in% "BP" ~ "#99d8c9",
    GO.cat %in% "MF" ~ "#3690c0"
  )) %>% mutate(GO.Rep = str_to_sentence(GO.Rep)) %>% 
  mutate(GO.Desc = str_to_sentence(GO.Desc)) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "dna", "DNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "ii", "II")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "ii", "II")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mapk", "MAPK")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mapk", "MAPK")) %>%
  print(width = Inf) 


#12 h


#Upregualtion
Gill_TFM_12_UP_GO_05_ALL_summary.1 <- Gill_TFM_12_UP_GO_05_ALL_summary %>% 
  filter(!is.na(GO.Rep)) %>%  # remove terms eliminated by REVIGO
  mutate(GO.cat = factor(GO.cat, levels = c("BP", "MF")),
         GO.Rep = factor(GO.Rep, levels = c('cellular response to transforming growth factor beta stimulus',
                                            'endosomal transport',
                                            'regulation of transcription from RNA polymerase II promoter',
                                            'cytokine receptor activity',
                                            'protein kinase activity',
                                            'protein tyrosine phosphatase activity',
                                            'RNA polymerase II regulatory region sequence-specific DNA binding',
                                            'transcription factor activity, RNA polymerase II core promoter proximal region sequence-specific binding',
                                            'transforming growth factor beta receptor binding'
         )),
         GO.Desc = factor(GO.Desc)) %>% 
  arrange(GO.cat, GO.Rep, desc(Gene.no), GO.Desc) %>%
  select( GO.cat, GO.Rep, GO.Desc, GO.Term, Gene.no, Genes) %>% 
  mutate(Colour = case_when(
    GO.cat %in% "BP" ~ "#99d8c9",
    GO.cat %in% "MF" ~ "#3690c0"
  )) %>% mutate(GO.Rep = str_to_sentence(GO.Rep)) %>% 
  mutate(GO.Desc = str_to_sentence(GO.Desc)) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "dna", "DNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "ii", "II")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "ii", "II")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mapk", "MAPK")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mapk", "MAPK")) %>%
  print(width = Inf) 

#Downregualtion
Gill_TFM_12_DOWN_GO_05_ALL_summary.1 <- Gill_TFM_12_DOWN_GO_05_ALL_summary %>% 
  filter(!is.na(GO.Rep)) %>%  # remove terms eliminated by REVIGO
  mutate(GO.cat = factor(GO.cat, levels = c("BP", "MF")),
         GO.Rep = factor(GO.Rep, levels = c('L-serine metabolism',
                                            'positive regulation of intracellular signal transduction'
         )),
         GO.Desc = factor(GO.Desc)) %>% 
  arrange(GO.cat, GO.Rep, desc(Gene.no), GO.Desc) %>%
  select( GO.cat, GO.Rep, GO.Desc, GO.Term, Gene.no, Genes) %>% 
  mutate(Colour = case_when(
    GO.cat %in% "BP" ~ "#99d8c9",
    GO.cat %in% "MF" ~ "#3690c0"
  )) %>% mutate(GO.Rep = str_to_sentence(GO.Rep)) %>% 
  mutate(GO.Desc = str_to_sentence(GO.Desc)) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "dna", "DNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "ii", "II")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "ii", "II")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mapk", "MAPK")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mapk", "MAPK")) %>%
  mutate(GO.Desc = str_replace(GO.Desc, "t cell", "T cell")) %>%
  mutate(GO.Desc = str_replace(GO.Desc, "FaT cell", "Fat cell")) %>%
  print(width = Inf) 



#24 


#Upregualtion
Gill_TFM_24_UP_GO_05_ALL_summary.1 <- Gill_TFM_24_UP_GO_05_ALL_summary %>% 
  filter(!is.na(GO.Rep)) %>%  # remove terms eliminated by REVIGO
  mutate(GO.cat = factor(GO.cat, levels = c("BP", "MF")),
         GO.Rep = factor(GO.Rep, levels = c('endosomal transport',
                                            'regulation of transcription from RNA polymerase II promoter',
                                            'venous blood vessel development',
                                            '2 iron, 2 sulfur cluster binding',
                                            'phosphatidylinositol kinase activity',
                                            'protein tyrosine phosphatase activity',
                                            'transcription regulatory region sequence-specific DNA binding',
                                            'transcriptional activator activity, RNA polymerase II transcription regulatory region sequence-specific binding',
                                            'transforming growth factor beta binding')),
         GO.Desc = factor(GO.Desc)) %>% 
  arrange(GO.cat, GO.Rep, desc(Gene.no), GO.Desc) %>%
  select( GO.cat, GO.Rep, GO.Desc, GO.Term, Gene.no, Genes) %>% 
  mutate(Colour = case_when(
    GO.cat %in% "BP" ~ "#99d8c9",
    GO.cat %in% "MF" ~ "#3690c0"
  )) %>% mutate(GO.Rep = str_to_sentence(GO.Rep)) %>% 
  mutate(GO.Desc = str_to_sentence(GO.Desc)) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "dna", "DNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "ii", "II")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "ii", "II")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mapk", "MAPK")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mapk", "MAPK")) %>%
  print(width = Inf) 



#Downregualtion
Gill_TFM_24_DOWN_GO_05_ALL_summary.1 <- Gill_TFM_24_DOWN_GO_05_ALL_summary %>% 
  filter(!is.na(GO.Rep)) %>%  # remove terms eliminated by REVIGO
  mutate(GO.cat = factor(GO.cat, levels = c("BP", "MF")),
         GO.Rep = factor(GO.Rep, levels = c('antigen processing and presentation of exogenous peptide antigen via MHC class I',
                                            'cytokine-mediated signaling pathway',
                                            'mRNA splicing, via spliceosome',
                                            'protein complex assembly',
                                            'SCF-dependent proteasomal ubiquitin-dependent protein catabolism',
                                            'viral life cycle',
                                            'cadherin binding',
                                            'G-rich strand telomeric DNA binding',
                                            'GTP binding',
                                            'peptidase activator activity',
                                            'peptide disulfide oxidoreductase activity',
                                            'proteasome-activating ATPase activity',
                                            'protein disulfide isomerase activity',
                                            'RNA-directed DNA polymerase activity',
                                            'small protein activating enzyme activity',
                                            'tumor necrosis factor-activated receptor activity'
         )),
         GO.Desc = factor(GO.Desc)) %>% 
  arrange(GO.cat, GO.Rep, desc(Gene.no), GO.Desc) %>%
  select( GO.cat, GO.Rep, GO.Desc, GO.Term, Gene.no, Genes) %>% 
  mutate(Colour = case_when(
    GO.cat %in% "BP" ~ "#99d8c9",
    GO.cat %in% "MF" ~ "#3690c0"
  )) %>% mutate(GO.Rep = str_to_sentence(GO.Rep)) %>% 
  mutate(GO.Desc = str_to_sentence(GO.Desc)) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Rna", "RNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "dna", "DNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "ii", "II")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtpase", "GTPase")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "MRNA", "mRNA")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "mhc", "MHC")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Mhc", "MHC")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "gtp", "GTP")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Gtp", "GTP")) %>%
  mutate(GO.Rep = str_replace(GO.Rep, "atp", "ATP")) %>% 
  mutate(GO.Rep = str_replace(GO.Rep, "Atp", "ATP")) %>%
  mutate(GO.Rep = str_replace(GO.Rep, "class i", "class I")) %>%
  mutate(GO.Rep = str_replace(GO.Rep, "class i", "class I")) %>%
  mutate(GO.Desc = str_replace(GO.Desc, "rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Rna", "RNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Dna", "DNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "ii", "II")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtpase", "GTPase")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "MRNA", "mRNA")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mhc", "MHC")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Mapk", "MAPK")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "mapk", "MAPK")) %>%
  mutate(GO.Desc = str_replace(GO.Desc, "class i", "class I")) %>%
  mutate(GO.Desc = str_replace(GO.Desc, "class i", "class I")) %>%
  mutate(GO.Desc = str_replace(GO.Desc, "gtp", "GTP")) %>% 
  mutate(GO.Desc = str_replace(GO.Desc, "Gtp", "GTP")) %>% 
  print(width = Inf) 


#Making the plots
#For a lot of these plots, I was filtering out specific GO terms on interest which involved
#the filter steps that are below. In most cases, MF was redundant and was put into the supplementary materials
#so that you will see a 'BP' vs 'MF' filtering upfront in most of the figures. 
#On that note, figures will be sperated between the manuscript (MS) or supplimentary figures


#To set an object for times new roman font for the figure text
windowsFonts(Times=windowsFont("TT Times New Roman")) 


#6 h


#UP

#Separating out BP from MF 
Gill6_UP_BP<- Gill_TFM_6_UP_GO_05_ALL_summary.1 %>% filter(GO.cat == "BP") %>% print()
Gill6_UP_MF<- Gill_TFM_6_UP_GO_05_ALL_summary.1 %>% filter(GO.cat == "MF") %>% print()


#getting my GO.Rep desired terms 

#Filtering out GO.Rep for MS figure
MS_fig_6U <- Gill6_UP_BP %>% 
  filter(GO.Rep =="Regulation of transcription from RNA polymerase II promoter") %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()



#Setting whats left for the supplementary figure

Joined_supp_6h_UP <- Gill_TFM_6_UP_GO_05_ALL_summary.1 %>% 
  filter(GO.Rep !="Regulation of transcription from RNA polymerase II promoter") %>% 
  arrange(GO.Rep, Gene.no) %>%  select(GO.Rep, GO.Desc, Gene.no, GO.cat, Colour) %>% 
  print()

Joined_supp_6h_UP.1<-Joined_supp_6h_UP %>% arrange( GO.cat) %>%  print()


Joined_supp_6h_UP.2$Colour<-as.factor(Joined_supp_6h_UP.2$Colour)


Joined_supp_6h_UP.2$Colour <- factor(Joined_supp_6h_UP.2$Colour)
Joined_supp_6h_UP.2$GO.Desc <- factor(Joined_supp_6h_UP.2$GO.Desc)
Joined_supp_6h_UP.2$GO.cat <- factor(Joined_supp_6h_UP.2$GO.cat)
Joined_supp_6h_UP.2$GO.Rep <- factor(Joined_supp_6h_UP.2$GO.Rep)
Joined_supp_6h_UP.2$Order <- factor(Joined_supp_6h_UP.2$Order)

str(Joined_supp_6h_UP.2)


Joined_supp_6h_UP.2$GO.Rep<-factor(Joined_supp_6h_UP.2$GO.Rep, levels = c('Cell junction assembly',
                                                                          'Heart development',
                                                                          'Regulation of transcription from RNA polymerase II promoter',
                                                                          'Chemokine receptor activity',
                                                                          'Phosphatidylinositol 3-kinase activity',
                                                                          'RNA polymerase II transcription factor binding',
                                                                          'Transcription factor activity, RNA polymerase II distal enhancer sequence-specific binding',
                                                                          'Transcription regulatory region DNA binding'))



#Manuscript figure

Gill_TFM_6_UP_GO_05_key_terms<- ggplot(MS_fig_6U)+
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill = "grey86")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title = element_blank(), 
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+
  scale_x_continuous(limits = c(0,70), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("6h")+
  theme(plot.title = element_text(hjust = 0.5, size = 21, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_6_UP_05_GOlist_BP ONLY_KEY fig_Sept_21.tiff",
         height =  6, width = 12)






#Supplementary figure
Gill_TFM_24_UP_supp_MF <- ggplot(Joined_supp_6h_UP.2) +
  geom_col(aes(x = Gene.no, y =reorder(GO.Desc, Gene.no), fill = Colour))+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19, family = "Times"),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,30), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  theme(legend.position = "None" )+
  scale_fill_manual(values = c("#99d8c9", "#3690c0"))+
  ggtitle("6h bluegill gill, upregulated")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_6_UP_05_GOlist_suppl_Sept_21.tiff",
         height = 8.5, width = 12)








#6 h down


Gill6_DOWN_BP<- Gill_TFM_6_DOWN_GO_05_ALL_summary.1 %>% filter(GO.cat == "BP") %>% print()

#MS filter
MS_fig_6D <- Gill6_DOWN_BP %>% 
  filter(GO.Rep =="Negative regulation of multicellular organismal process") %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()



#Suppl filter
Suppl_fig_6D <- Gill_TFM_6_DOWN_GO_05_ALL_summary.1 %>% 
  filter(GO.Rep !="Negative regulation of multicellular organismal process") %>% 
  arrange(GO.Rep, Gene.no) %>%  select(GO.Rep, GO.Desc, Gene.no, GO.cat, Colour) %>% 
  print()

Suppl_fig_6D$Colour <- factor(Suppl_fig_6D$Colour)
Suppl_fig_6D$GO.Desc <- factor(Suppl_fig_6D$GO.Desc)
Suppl_fig_6D$GO.cat <- factor(Suppl_fig_6D$GO.cat)
Suppl_fig_6D$GO.Rep <- factor(Suppl_fig_6D$GO.Rep)

Suppl_fig_6D$GO.Rep<-factor(Suppl_fig_6D$GO.Rep, levels = c('Actin cytoskeleton reorganization',
                                                            'Negative regulation of multicellular organismal process',
                                                            'Proteoglycan biosynthesis',
                                                            'G-protein coupled receptor binding',
                                                            'Transmembrane receptor protein tyrosine kinase activity'))

Suppl_fig_6D.1 <-Suppl_fig_6D %>% arrange(GO.cat) %>% mutate(erp = Colour) %>% print()

#MS plot 

Gill_TFM_6_DOWN_GO_05_key_terms<- ggplot(MS_fig_6D)+
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill = "grey86")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title = element_blank(), 
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+
  scale_x_continuous(limits = c(0,20), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("6h")+
  theme(plot.title = element_text(hjust = 0.5, size = 21, face = "bold", family = "Times"))+
  ggsave("/BG_Gill_TFM_6_DOWN_05_GOlist_BP ONLY_KEY fig_Sept_21.tiff",
         height =  4, width = 12)



#Suppl plot

Gill_TFM_24_UP_supp_MF <- ggplot(Suppl_fig_6D.1) +
  geom_col(aes(x = Gene.no, y =reorder(GO.Desc, Gene.no), fill = erp))+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19, family = "Times"),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,30), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  theme(legend.position = "None" )+
  scale_fill_manual(values = c('#3690c0','#99d8c9' ))+
  ggtitle("6h bluegill gill, downregulated")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_6_DOWN_05_GOlist_suppl_Sept_21.tiff",
         height =  4, width = 12)




#12 h



#12 UP

#splitting the BP and MP terms up for ease

Gill12_UP_BP<- Gill_TFM_12_UP_GO_05_ALL_summary.1 %>% filter(GO.cat == "BP") %>% print()
Gill12_UP_MF<- Gill_TFM_12_UP_GO_05_ALL_summary.1 %>% filter(GO.cat == "MF") %>% print()


#MS filter
MS_fig_12U <- Gill12_UP_BP %>% 
  filter(GO.Rep == "Cellular response to transforming growth factor beta stimulus"|GO.Rep == "Regulation of transcription from RNA polymerase II promoter") %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()


#Suppl filter
Suppl_fig_12U <- Gill_TFM_12_UP_GO_05_ALL_summary.1 %>% 
  filter(GO.Rep != "Cellular response to transforming growth factor beta stimulus") %>%
  filter(GO.Rep != "Regulation of transcription from RNA polymerase II promoter") %>% 
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no, GO.cat, Colour) %>% 
  mutate(erp = Colour) %>% select(-Colour) %>% 
  print()


Suppl_fig_12U$erp <- factor(Suppl_fig_12U$erp)
Suppl_fig_12U$GO.Desc <- factor(Suppl_fig_12U$GO.Desc)
Suppl_fig_12U$GO.cat <- factor(Suppl_fig_12U$GO.cat)
Suppl_fig_12U$GO.Rep <- factor(Suppl_fig_12U$GO.Rep)

Suppl_fig_12U$GO.Rep<-factor(Suppl_fig_12U$GO.Rep, levels = c(
  'Endosomal transport',
  
  'Cytokine receptor activity',
  'Protein kinase activity',
  'Protein tyrosine phosphatase activity',
  'RNA polymerase II regulatory region sequence-specific DNA binding',
  'Transcription factor activity, RNA polymerase II core promoter proximal region sequence-specific binding',
  'Transforming growth factor beta receptor binding'))

Suppl_fig_12U.1 <-Suppl_fig_12U %>% arrange(GO.cat)  %>% print()


#MS plot 
Gill_TFM_12_UP_GO_05_key_terms<- ggplot(MS_fig_12U)+
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill = "grey34")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(italic("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title = element_blank(),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+
  scale_x_continuous(limits = c(0,150), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("12h")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_12_UP_05_GOlist_BP ONLY_KEY fig_Sept_21.tiff",
         height =  8.5, width = 12)


#Suppl plot
ggplot(Suppl_fig_12U.1) +
  geom_col(aes(x = Gene.no, y =reorder(GO.Desc, Gene.no), fill = Suppl_fig_12U.1$erp))+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19, family = "Times"),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,60), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  theme(legend.position = "None" )+
  scale_fill_manual(values = c('#3690c0','#99d8c9' ))+
  ggtitle("12h bluegill gill, upregulated")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_12_UP_05_GOlist_suppl_Sept_21.tiff",
         height = 11, width = 12)





#12 h down

Gill12_DOWN_BP<- Gill_TFM_12_DOWN_GO_05_ALL_summary.1 %>% filter(GO.cat == "BP") %>% print()

#MS filter
MS_fig_12D <- Gill12_DOWN_BP %>% 
  filter(GO.Rep == "Positive regulation of intracellular signal transduction") %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()


#Suppl filter
Suppl_fig_12D <- Gill12_DOWN_BP %>% 
  filter(GO.Rep != "Positive regulation of intracellular signal transduction") %>% 
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()


#MS plot
Gill_TFM_12_DOWN_GO_05_key_terms<- ggplot(MS_fig_12D)+
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill = "grey34")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title = element_blank(), 
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+
  scale_x_continuous(limits = c(0,60), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("12h")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_12_DOWN_05_GOlist_BP ONLY_KEY fig_Sept_21.tiff",
         height =  6, width = 12)



#Suppl plot
ggplot(Suppl_fig_12D) +
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill  = "#99d8c9")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19, family = "Times"),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+
  scale_x_continuous(limits = c(0,5), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("12h bluegill gill, downregulated, biological processes")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  scale_fill_manual(values = c("BP" = "#99d8c9",  "MF" = "#3690c0"), 
                    labels = c( "BP" = "Biological Processes", "MF" = "Molecular Function"))+
  ggsave("BG_Gill_TFM_12_DOWN_05_GOlist_suppl_Sept_21.tiff",
         height = 2.5, width = 12)





#24 h




#24 UP

Gill24_UP_BP<- Gill_TFM_24_UP_GO_05_ALL_summary.1 %>% filter(GO.cat == "BP") %>% print()
Gill24_UP_MF<- Gill_TFM_24_UP_GO_05_ALL_summary.1 %>% filter(GO.cat == "MF") %>% print()

#MS filter

MS_fig_24U <- Gill24_UP_BP %>% 
  filter(GO.Rep ==c("Regulation of transcription from RNA polymerase II promoter")) %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()

#Suppl filter

Suppl_fig_24U <- Gill_TFM_24_UP_GO_05_ALL_summary.1 %>% 
  filter(GO.Rep != "Regulation of transcription from RNA polymerase II promoter") %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no, GO.cat, Colour) %>% 
  mutate(erp = Colour) %>% select(-Colour) %>% 
  print()


Suppl_fig_24U$erp <- factor(Suppl_fig_24U$erp)
Suppl_fig_24U$GO.Desc <- factor(Suppl_fig_24U$GO.Desc)
Suppl_fig_24U$GO.cat <- factor(Suppl_fig_24U$GO.cat)
Suppl_fig_24U$GO.Rep <- factor(Suppl_fig_24U$GO.Rep)

Suppl_fig_24U$GO.Rep<-factor(Suppl_fig_24U$GO.Rep, levels = c(
  'Endosomal transport',
  
  'Venous blood vessel development',
  '2 iron, 2 sulfur cluster binding',
  'Phosphatidylinositol kinase activity',
  'Protein tyrosine phosphatase activity',
  'Transcription regulatory region sequence-specific DNA binding',
  'Transcriptional activator activity, RNA polymerase II transcription regulatory region sequence-specific binding',
  'Transforming growth factor beta binding'))

Suppl_fig_24U.1 <-Suppl_fig_24U %>% arrange(GO.cat)  %>% print()


#MS plot

Gill_TFM_24_UP_GO_05_key_terms<- ggplot(MS_fig_24U)+
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill = "grey2")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(italic("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19), 
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,180), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("24h")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_24_UP_05_GOlist_BP ONLY_KEY fig_Sept_21.tiff",
         height =  6, width = 12)



#Suppl plot

ggplot(Suppl_fig_24U.1) +
  geom_col(aes(x = Gene.no, y =reorder(GO.Desc, Gene.no), fill = Suppl_fig_24U.1$erp))+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19, family = "Times"),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,60), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  theme(legend.position = "None" )+
  scale_fill_manual(values = c('#3690c0','#99d8c9' ))+
  ggtitle("24h bluegill gill, upregulated")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_24_UP_05_GOlist_suppl_Sept_21.tiff",
         height = 10, width = 12)






#24 down



D24_TFM_BP <- Gill_TFM_24_DOWN_GO_05_ALL_summary.1 %>% filter(GO.cat == "BP") %>% print()
D24_TFM_MF <- Gill_TFM_24_DOWN_GO_05_ALL_summary.1 %>% filter(GO.cat == "MF") %>% print()


#MS filter
MS_fig_24D <- D24_TFM_BP %>% 
  filter(GO.Rep =='Antigen processing and presentation of exogenous peptide antigen via MHC class I'| GO.Rep == 
           'Cytokine-mediated signaling pathway'| GO.Rep =='Scf-dependent proteasomal ubiquitin-dependent protein catabolism'| GO.Rep ==
           'Viral life cycle') %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no) %>% 
  print()



#Suppl filter

Suppl_fig_24D<-Gill_TFM_24_DOWN_GO_05_ALL_summary.1 %>% 
  filter(GO.Rep !='Antigen processing and presentation of exogenous peptide antigen via MHC class I') %>% 
  filter(GO.Rep != 'Cytokine-mediated signaling pathway') %>% 
  filter(GO.Rep !='Scf-dependent proteasomal ubiquitin-dependent protein catabolism') %>% 
  filter(GO.Rep !=
           'Viral life cycle') %>%
  arrange(GO.Rep, Gene.no) %>% select(GO.Rep, GO.Desc, Gene.no, GO.cat, Colour) %>% 
  mutate(erp = Colour) %>% select(-Colour) %>% 
  print()

View(Gill_TFM_24_DOWN_GO_05_ALL_summary.1)
Suppl_fig_24D$erp <- factor(Suppl_fig_24D$erp)
Suppl_fig_24D$GO.Desc <- factor(Suppl_fig_24D$GO.Desc)
Suppl_fig_24D$GO.cat <- factor(Suppl_fig_24D$GO.cat)
Suppl_fig_24D$GO.Rep <- factor(Suppl_fig_24D$GO.Rep)

Suppl_fig_24D$GO.Rep<-factor(Suppl_fig_24D$GO.Rep, levels = c(
  'mRNA splicing, via spliceosome',
  'Protein complex assembly',
  'Scf-dependent proteasomal ubiquitin-dependent protein catabolism',
  'viral life cycle',
  'Cadherin binding',
  'G-rich strand telomeric DNA binding',
  'GTP binding',
  'Peptidase activator activity',
  'Peptide disulfide oxidoreductase activity',
  'Proteasome-activating ATPase activity',
  'Protein disulfide isomerase activity',
  'RNA-directed DNA polymerase activity',
  'Small protein activating enzyme activity',
  'Tumor necrosis factor-activated receptor activity'))
Suppl_fig_24D
Suppl_fig_24D.1 <-Suppl_fig_24D %>% arrange(GO.cat)  %>% print()


#MS plot

Gill_TFM_24_DOWN_GO_05_key_terms<- ggplot(MS_fig_24D)+
  geom_col(aes(x = Gene.no, y = reorder(GO.Desc, Gene.no)), fill = "grey2")+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(italic("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19), 
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,180), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  ggtitle("24h")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_24_DOWN_05_GOlist_BP ONLY_KEY fig_Sept_21.tiff",
         height =  13, width = 12)



#Suppl plot 

Gill_TFM_24_D_GO_05_suppl <- ggplot(Suppl_fig_24D.1) +
  geom_col(aes(x = Gene.no, y =reorder(GO.Desc, Gene.no), fill = Suppl_fig_24D.1$erp))+
  facet_grid(rows = vars(GO.Rep), scales = "free", space = "free", 
             labeller = label_wrap_gen(50))+
  scale_y_discrete(labels = label_wrap_gen(40))+
  labs(x = expression(bold("Gene number")),
       y = element_blank())+
  theme(panel.background = element_blank(),
        strip.background = element_rect(fill = "0.5"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank()) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 13, family = "Times"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 19, family = "Times"),  
        axis.text.x = element_text(size = 21, family = "Times"),
        axis.text.y = element_text(size = 10, family = "Times"))+
  scale_x_continuous(limits = c(0,200), expand = c(0, 0))+
  theme(axis.ticks.y = element_line(), axis.ticks.length = unit(5, "pt"))+
  theme(plot.margin=unit(c(.6,.6,.6,.6),"cm"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  theme(legend.position = "None" )+
  scale_fill_manual(values = c('#3690c0','#99d8c9' ))+
  ggtitle("24h bluegill gill, downregulated")+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Times"))+
  ggsave("BG_Gill_TFM_24_DOWN_05_GOlist_suppl_Sept_21.tiff",
         height = 16, width = 12)




#####Heat maps #####

library(tidyverse)
library(superheat)

setwd("")

#Making heatmaps for both treatment level effects and individual fish expression patterns
#Used the R package 'superheat' to make them

#Treatment level heatmaps


#Read in the annotated results file
all.results_annot<-read.csv("all.results_annot.csv")


dim(all.results_annot) #46163 56 
names(all.results_annot)

#Making a simplified tibble by extracting out the clusters that showed significant differential expression
#at one of the timepoints 


heat_data_TFM<- all.results_annot %>%  
  select(Cluster, sp.BX_gene.name, C.6vT.6_logFC,C.6vT.6_FDR, C.12vT.12_logFC,C.12vT.12_FDR, 
         C.24vT.24_logFC,C.24vT.24_FDR) %>% filter (C.6vT.6_FDR <0.05 | C.12vT.12_FDR <0.05 | C.24vT.24_FDR <0.05) %>% 
  rename("6h" = C.6vT.6_logFC, "12h" =C.12vT.12_logFC, "24h" = C.24vT.24_logFC)
dim(heat_data_TFM) #5810 8

#Getting cluster names associated with a row
heat_data_TFM_ALL <- column_to_rownames(as.data.frame(heat_data_TFM), var = "Cluster") %>% 
  select("6h", "12h", "24h")
head(heat_data_TFM_ALL)

#making the figure

jpeg("heatmap_BG_TFM_gills_Overall.jpeg",width=3.5, height=5.5, units="in", res=300, quality=100)

superheat(heat_data_TFM_ALL,
          pretty.order.rows = TRUE,
          grid.hline.col = "white",
          grid.vline.col = "white",
          bottom.label.text.size = 3,
          bottom.label.col = "white",
          left.label.text.size = 3,
          left.label.text.alignment = "center",
          left.label.col = "black",
          legend.width = 1,
          legend.text.size = 10,
          padding = 0.25,
          scale = T
          )


dev.off()



#Individual heatmaps


#Read in extracted Log CPM values

Log_CPM_TFM <- read.csv("Log CPM_norm values_BG_TFM_gill.csv")


#Filtering out genes that were signification the treatment heat map that are in the individual fish


sig_clusters <- heat_data_TFM %>%  select(Cluster) %>% print()
dim (sig_clusters) #5810 1


Log_CPM_TFM.1 <- Log_CPM_TFM %>% filter(Cluster %in% sig_clusters$Cluster) %>% print()
dim(Log_CPM_TFM.1) #5810 40

head(Log_CPM_TFM)


# X             Cluster      C6.1G        C6.4G      C6.5G     C6.7G    C12.8G        T6.1G     T12.3G      T24.5G      C6.8G     C24.1G      T6.2G     T12.4G     T24.7G
# 1 Cluster_51425.34829 Cluster_51425.34829 1.86622209  1.646535230  1.7386442 1.5597811 1.4756243  1.122600666  0.7596749  1.76211741  0.9721919  1.4672314  1.6287417  1.5764375  0.4517979
# 2     Cluster_37848.2     Cluster_37848.2 4.87162377  5.017308148  4.8820555 4.9160345 4.6089723  4.822955726  4.9867003  4.82477289  4.8396613  4.4272476  4.7276579  4.8218727  4.6893232
# 3 Cluster_51425.14595 Cluster_51425.14595 0.48072474  0.001246198 -0.4133579 0.5498922 0.6918796 -0.409239459 -0.6693054 -0.03576739 -0.2837776 -3.0593805 -0.4363974 -1.5400478 -0.4850683
# 4  Cluster_51425.9311  Cluster_51425.9311 1.10580537  1.587613460  0.9262714 2.1213088 0.3215722 -0.004047926  1.4817563  0.76036128  0.7726465  1.6673673  1.2026817  2.5615625  1.3954755
# 5 Cluster_51425.23955 Cluster_51425.23955 0.05239663 -0.414494195 -1.1093510 0.5697892 0.1538781 -0.489394655 -1.3941348 -0.03576739 -0.6634339 -0.3836377 -0.5118504  0.2533208 -0.3997011
# 6      Cluster_7525.2      Cluster_7525.2 5.11571277  5.096048279  5.2301408 5.0018386 4.9809175  5.009493706  4.9407781  5.08284679  4.9765942  4.8063283  4.8513722  4.9372947  5.0664827


head(Log_CPM_TFM.1)


# X             Cluster    C6.1G    C6.4G    C6.5G    C6.7G   C12.8G    T6.1G   T12.3G   T24.5G    C6.8G   C24.1G    T6.2G   T12.4G   T24.7G   C12.1G   C24.2G    T6.4G   T12.5G
# 1     Cluster_37848.2     Cluster_37848.2 4.871624 5.017308 4.882056 4.916035 4.608972 4.822956 4.986700 4.824773 4.839661 4.427248 4.727658 4.821873 4.689323 4.756112 4.676073 4.796091 4.908505
# 2 Cluster_51425.11981 Cluster_51425.11981 5.162990 5.165685 4.870764 5.122096 4.825625 5.980493 6.089697 5.781000 4.947268 4.051147 5.409985 5.646667 5.529283 4.929524 5.273047 5.888954 5.749825
# 3 Cluster_51425.32231 Cluster_51425.32231 3.853922 4.014450 4.029134 3.926691 3.944626 3.874911 4.127651 4.023068 3.767555 3.859993 3.858178 4.136522 4.231782 4.059065 3.734371 3.995559 4.183232
# 4    Cluster_31777.15    Cluster_31777.15 5.124799 5.212749 5.262903 5.411141 5.328167 5.108117 5.679593 5.420845 5.332397 5.199193 5.245738 5.639867 5.214279 5.026173 4.987156 5.282707 5.651826
# 5     Cluster_39959.2     Cluster_39959.2 5.681470 5.732934 5.726217 5.592280 5.700159 5.869388 5.702379 5.741040 5.670582 5.705869 5.750921 5.648023 5.660539 5.628701 5.479013 5.704085 5.748655
# 6     Cluster_35915.1     Cluster_35915.1 3.859870 4.042144 4.192235 3.909227 3.653775 3.378241 3.534945 3.844377 4.214180 4.794155 3.534023 3.748665 3.898443 3.879667 4.376363 3.252530 3.414512

#removing the extra column with the cluster ID's
Log_CPM_TFM.1 <- Log_CPM_TFM.1 %>% select(-X) %>% print()





heat_data_TFM_indy<- Log_CPM_TFM.1  


dim(heat_data_TFM_indy) #5810 39


Indy_heat_data_TFM_ALL <- column_to_rownames(as.data.frame(heat_data_TFM_indy), var = "Cluster") %>% 
  select(everything())
head(Indy_heat_data_TFM_ALL)

#Extracting the column names to make better labels for plot

Colnames<-colnames(Indy_heat_data_TFM_ALL)
write.csv(Colnames, "Heatmap_TFM_indy_colnames.csv")

#Renaming the fish names here for clairty on the figure

Indy_heat_data_TFM_ALL.1 <- Indy_heat_data_TFM_ALL %>% rename('Control 6h.1' = 'C6.1G',
                                                              'Control 6h.4' = 'C6.4G',
                                                              'Control 6h.5' = 'C6.5G',
                                                              'Control 6h.7' = 'C6.7G',
                                                              'Control 12h.8' = 'C12.8G',
                                                              'TFM 6h.1' = 'T6.1G',
                                                              'TFM 12h.3' = 'T12.3G',
                                                              'TFM 24h.5' = 'T24.5G',
                                                              'Control 6h.8' = 'C6.8G',
                                                              'Control 24h.1' = 'C24.1G',
                                                              'TFM 6h.2' = 'T6.2G',
                                                              'TFM 12h.4' = 'T12.4G',
                                                              'TFM 24h.7' = 'T24.7G',
                                                              'Control 12h.1' = 'C12.1G',
                                                              'Control 24h.2' = 'C24.2G',
                                                              'TFM 6h.4' = 'T6.4G',
                                                              'TFM 12h.5' = 'T12.5G',
                                                              'TFM 24h.8' = 'T24.8G',
                                                              'Control 12h.2' = 'C12.2G',
                                                              'Control 24h.4' = 'C24.4G',
                                                              'TFM 6h.6' = 'T6.6G',
                                                              'TFM 12h.6' = 'T12.6G',
                                                              'TFM 24h.9' = 'T24.9G',
                                                              'Control 12h.3' = 'C12.3G',
                                                              'TFM 6h.7' = 'T6.7G',
                                                              'TFM 12h.7' = 'T12.7G',
                                                              'Control 12h.4' = 'C12.4G',
                                                              'Control 24h.6' = 'C24.6G',
                                                              'TFM 6h.9' = 'T6.9G',
                                                              'TFM 12h.8' = 'T12.8G',
                                                              'Control 12h.5' = 'C12.5G',
                                                              'Control 24h.7' = 'C24.7G',
                                                              'TFM 12h.1' = 'T12.1G',
                                                              'TFM 12h.9' = 'T12.9G',
                                                              'Control 12h.6' = 'C12.6G',
                                                              'Control 24h.8' = 'C24.8G',
                                                              'TFM 12h.2' = 'T12.2G',
                                                              'TFM 24h.4' = 'T24.4G') %>% print()



#Reading in an rrdering vector to reorganize the cols in the plot (was done in excel)

Order_vector <- read.csv("heatmap_matrix_TFM_gill_BG.csv", header=T) %>% print()

dim(Order_vector) #38

OG_order <- colnames(Indy_heat_data_TFM_ALL.1)
OG_order.1 <- OG_order %>% as.tibble %>% rename("fish" = "value") %>% print()



Order_vector.1 <- left_join(OG_order.1, Order_vector, by = "fish") %>% print()




#Making the individual heatmap

jpeg("heatmap_BG_TFM_gill_individualfish_simple.jpeg",width=10, height=10, units="in", res=300, quality=100)

superheat(Indy_heat_data_TFM_ALL.1,
          pretty.order.rows = TRUE,
          grid.hline.col = "white",
          grid.vline.col = "white",
          bottom.label.text.size = 4,
          bottom.label.col = "white",
          left.label.text.size = 3,
          left.label.text.alignment = "center",
          left.label.col = "black",
          legend.width = 1,
          legend.text.size = 9,
          padding = 0.25,
          scale = T,
          bottom.label.text.angle = 90,
          order.cols = order(Order_vector.1$order), #sorts the col's by my provided vector
          column.title = "Fish ID",
          column.title.size = 6) 


dev.off()



######Detoxification Gene Filtering#######


#As noted in the manuscript, we were interested in getting genes involved in TFM detoxification for Phase I-III 
#biotransformation processes
#The below outlines the filtering processes to get the raw gene names/id's which were then manually checked for their
#function using Uniprot descriptions to do so (https://www.uniprot.org/)


library(tidyr)
library(tidyverse)
library(car)



####UGT isolation#####

#Pulling out any of the UGT's in the transcriptomes for presence/absence independent of DE patterns
#Just inherent UGT gene availability 

#read in the transcriptome and make all genes in upper case 

BGtranscriptome <- read_csv("bluegill_transcriptome.csv") %>% mutate(gene.name = toupper(gene.name)) %>% print()

#filter out all UGT's in transcriptome and add column with species ID
UGT_BG <- dplyr::filter(BGtranscriptome, grepl('UGT', gene.name)) %>% add_column(Exp.Species = 'BG') %>% 
  select(Exp.Species, everything()) %>% print() #19 hits



#Writing out file and will curate in Excel/Uniprot
UGT_BG %>% write.csv("BG_Gill_UGTisoforms transcriptome June 24.csv")


#Getting those UGT genes out to see if any had DE in the gills
#Using the compiled Cluster files from the gene counts in each of the files and making the gene names in upper case

BG_gill <- read_csv("All_clusters_significant_TFM_gills.csv") %>% select(-X1, -Treatment) %>% mutate(sp.BX_gene.name = toupper(sp.BX_gene.name)) %>% print()




#Filtering out differently expressed UGT genes in the bluegill gill

UGT_BG_gill <-dplyr::filter(BG_gill, grepl('UGT', sp.BX_gene.name)) %>% add_column(Exp.Species = 'BG') %>%  
  add_column(Tissue = 'Gill') %>%  select(Exp.Species,Tissue, everything())%>% print() 

#3 hits

#Write out the file 
UGT_BG_gill %>% write.csv("BG_Gill_Gill_UGTs_Sig_June 24.csv")



#####Phase I genes isolation#####


Phase1_BG_gill <- dplyr::filter(BG_gill, grepl('CYP|ADH|MAO|PON|ALDH', sp.BX_gene.name)) %>% add_column(Exp.Species = 'BG') %>%  #
  add_column(Tissue = 'Gill') %>%  select(Exp.Species,Tissue, everything())%>% print() 
#23 genes



#Write out the file 
Phase1_BG_gill %>% write.csv("BG_Gill_Phase_I_sig_genes_June 24.csv")



#####Phase II enzymes ioslation#####



Phase2_BG_gill <- dplyr::filter(BG_gill, grepl('SULT|GST|GLYAT|NAT|MT', sp.BX_gene.name)) %>% add_column(Exp.Species = 'BG') %>%  #
  add_column(Tissue = 'Gill') %>%  select(Exp.Species,Tissue, everything())%>% print() 
#75 genes


Phase2_BG_gill %>% write.csv("BG_Gill_Phase_II_sig_genes_June 24.csv")



##### Phase III enzymes#####


Phase3_BG_gill <- dplyr::filter(BG_gill, grepl('SLC|ABC', sp.BX_gene.name)) %>% add_column(Exp.Species = 'BG') %>%  #
  add_column(Tissue = 'Gill') %>%  select(Exp.Species,Tissue, everything())%>% print() 
#167 genes


Phase3_BG_gill %>% write.csv("BG_Gill_Phase_III_sig_genes_June 21.csv")



























