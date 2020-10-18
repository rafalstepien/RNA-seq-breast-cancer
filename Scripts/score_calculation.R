library(glue)
homepath <- Sys.getenv('HOME')

cancer_kallisto_counts <- read.csv(glue('{homepath}/RNA-seq/RNA-seq-breast-cancer/Outputs/concatenate_files_output/cancer_kallisto_counts.csv'))
normal_kallisto_counts <- read.csv(glue('{homepath}/RNA-seq/RNA-seq-breast-cancer/Outputs/concatenate_files_output/normal_kallisto_counts.csv'))
cancer_kallisto_counts = cancer_kallisto_counts[, c(-6, -7)]

library(biomaRt)
mart <- useMart(biomart = 'ensembl', dataset='hsapiens_gene_ensembl')
cancer_ids = cancer_kallisto_counts$target_id

res_cancer <- getBM(attributes= c('ensembl_transcript_id', 'ensembl_gene_id', 'external_transcript_name', 'external_gene_name'),
             filters = 'ensembl_transcript_id',
             values = cancer_ids,
             mart=mart
             )
colnames(cancer_kallisto_counts)[1] = "ensembl_transcript_id"
cancer_counts = merge(cancer_kallisto_counts, res_cancer, by="ensembl_transcript_id", all.x = TRUE)
cancer_counts = cancer_counts[,c(1, 6, 7, 8, 2, 3, 4, 5)]

library(tidyverse)
reference_group_cancer = cancer_counts %>% filter(external_gene_name == 'ACTB' | 
                                            external_gene_name == 'GAPDH' |
                                            external_gene_name == 'RPLPO' |
                                            external_gene_name == 'GUS' |
                                            external_gene_name == 'TFRC')

GRB7_group_cancer =  cancer_counts %>% filter(external_gene_name == 'GRB7' | 
                                              external_gene_name == 'HER2')

ER_group_cancer = cancer_counts %>% filter(external_gene_name == 'ER' | 
                                           external_gene_name == 'PGR' |
                                           external_gene_name == 'BCL' |
                                           external_gene_name == 'SCUBE2')

proliferation_cancer =  cancer_counts %>% filter(external_gene_name == 'Ki67' | 
                                                  external_gene_name == 'STK15' |
                                                  external_gene_name == 'Survivin' |
                                                  external_gene_name == 'CCNB1' |
                                                  external_gene_name == 'MYBL2')

invasion_cancer = cancer_counts %>% filter(external_gene_name == 'MMP11' | 
                                           external_gene_name == 'CTSL2')

normal_ids = normal_kallisto_counts$target_id
res_normal <- getBM(attributes= c('ensembl_transcript_id', 'ensembl_gene_id', 'external_transcript_name', 'external_gene_name'),
                    filters = 'ensembl_transcript_id',
                    values = normal_ids,
                    mart=mart
)
colnames(normal_kallisto_counts)[1] = "ensembl_transcript_id"
normal_counts = merge(normal_kallisto_counts, res_normal, by="ensembl_transcript_id", all.x = TRUE)
normal_counts = normal_counts[,c(1, 6, 7, 8, 2, 3, 4, 5)]

reference_group_normal = normal_counts %>% filter(external_gene_name == 'ACTB' | 
                                                  external_gene_name == 'GAPDH' |
                                                  external_gene_name == 'RPLPO' |
                                                  external_gene_name == 'GUS' |
                                                  external_gene_name == 'TFRC')

GRB7_group_normal =  normal_counts %>% filter(external_gene_name == 'GRB7' | 
                                              external_gene_name == 'HER2')

ER_group_normal = normal_counts %>% filter(external_gene_name == 'ER' | 
                                           external_gene_name == 'PGR' |
                                           external_gene_name == 'BCL' |
                                           external_gene_name == 'SCUBE2')

proliferation_normal =  normal_counts %>% filter(external_gene_name == 'Ki67' | 
                                                  external_gene_name == 'STK15' |
                                                  external_gene_name == 'Survivin' |
                                                  external_gene_name == 'CCNB1' |
                                                  external_gene_name == 'MYBL2')

invasion_normal = normal_counts %>% filter(external_gene_name == 'MMP11' | 
                                           external_gene_name == 'CTSL2')


