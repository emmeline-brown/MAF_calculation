#Calculate minor allele frequency of a biallelic variant in Plink format
#Step 1: use Plink to select a single variant
#This is demonstrated on the rs53576 chr3:8804371-8804371 which  is a silent G to A change in the oxytocin receptor (OXTR)
#A is the minor allele

#load package

library(tidyverse)

#input files

my_var_map <- read_tsv("/path/to/snp/rs53576.map", col_names = F)
my_var_ped <- read_tsv("/path/to/snp/rs53576.ped", col_names = F)

#select relevant variable from .MAP file
my_var_map <-
  my_var_map %>%
  select(X2) %>% 
  rename(variant = X2)

#select relevant variable from .PED file
my_var_ped <-
  my_var_ped %>% 
  select(X7) %>% 
  rename(variant = X7)

#pass the rsID to the column name
map_ped <- my_var_map %>% 
  full_join(my_var_ped)
colnames(map_ped) <- map_ped[1, ] #this makes the header the first row
map_ped <- map_ped[-1, ] #this removes the first row

#remove surplus files
rm(my_var_ped, my_var_map)

#remove white spaces from variant
map_ped$rs53576 <- gsub(" ", "", map_ped$rs53576)

count_rs <- dplyr::summarise(group_by(map_ped, rs53576, count =n()))
count_rs <- 
  count_rs %>% 
  spread(key = rs53576, value=count)

#MAF calculations
#Search online for which allele is most common in the population, and which is the minor allele
MAF <- count_rs %>% 
  mutate(total_count = AA + GA + GG, #change these to other alleles if required
         minor_hom = AA *2) %>% #change the minor homozygote allele here
  mutate(allele_number = total_count * 2,
         allele_count = minor_hom + GA) %>% #change this to the heterozygote allele required
  mutate(MAF = allele_count/allele_number) %>% 
  mutate(total_individuals = AA + GA + AA) %>% #change if required
  select(-allele_number, -allele_count, -minor_hom, -total_count, -AA, -GA, -GG)

#write table
write_csv(MAF, "path/to/saving/SNP_MAF.csv")

#remove additional
rm(count_rs, MAF, map_ped)


