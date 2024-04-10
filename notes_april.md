FST:

```console
fst_main<-read.table("/project/ctb-grego/kplanedin/blue_deacon_widow/vcf/variant_only/filtered/filt0/fst/fst_sexed_deacon_only_allsamplesvcf_fst.txt",header=T)
fst_main<-fst_main %>%
  transform(chromosome=str_replace(chromosome,"genome_hic_scaffold_","")) %>%
  mutate(chromosome = as.numeric(chromosome)) %>% filter(chromosome < 25)
  
fst_filt<-fst_main %>% filter(!is.na(avg_wc_fst))

chr_lengths <- read_tsv("/project/ctb-grego/kplanedin/other_species/sebastes/Ssc.HIC.genome.fa.fai",
                        col_names = c("chromosome","length","bits","spacer1","spacer2"))

chr_lengths$chromosome <- sub("^genome_hic_scaffold_", "", chr_lengths$chromosome)
chr_lengths$chromosome <- as.numeric(chr_lengths$chromosome)
chr_lengths<-chr_lengths %>% filter(chromosome < 25)

chr_lengths %>%
  # Select only the columns we need
  dplyr::select(chromosome,length) %>%
  # Calculate cumulative position of each chromosome
  mutate(total=cumsum(length)-length) %>%
  # Remove the length column, we don't need it anymore.
  dplyr::select(-length) %>%
  # Add this info to the initial dataset
  left_join(fst_filt, ., by=c("chromosome"="chromosome")) %>%
  # Make sure that positions are still sorted
  arrange(chromosome, window_pos_1) %>%
  # Now create a new column with the cumulative position
  mutate( cumulative_pos=window_pos_1+total) -> fst_cumulative


plot_genome_fst<-fst_cumulative %>%
  ggplot(.) +
  geom_point( aes(x=cumulative_pos, y=avg_wc_fst,color=as.factor(chromosome)), 
              alpha=1, size=1) +
  # Alternate chromosome point colors
  scale_color_manual(values = rep(c("grey85", "grey"), nrow(chr_lengths) )) +
  # Custom X-axis
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center ) +
  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(title="Genome Wide FST of Male and Female Deacon (10 kbp)",x=NULL,y=NULL)
```

##################################
messy down here, but thats okay
##################################

gwas <- read_table("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/results/gwas/salmon.gwas.P1.assoc.linear")
gwas <- gwas %>%
  filter(TEST == "ADD")
gwas <- gwas %>%
  mutate(log_p = -log10(P))
gwas %>%
  filter(log_p > 1) %>%
  ggplot(.,aes(x=BP,y=log_p)) +
  geom_point()
gwas %>%
  arrange(desc(log_p)) %>%
  head()

file_list<-read.table("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/fst/combo_list.txt")
file_list <- unlist(file_list)
snp_list<-c(3185216,3185230,3185308,3185255)

snp_fst_list<-data.frame()


for (i in file_list){
  
combo_name_tmp<- gsub(paste0(".txt.weir.fst", "$"), "", i)
combo_name<- gsub(paste0("^", "fst_"), "", combo_name_tmp)

print(paste("Starting: ",combo_name))  

path<-paste0("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/fst/",i)


fst<-read.table(path,header=T)
fst_select<-fst %>% mutate(populations = combo_name) %>%
  filter(POS %in% snp_list)

snp_fst_list<-rbind(snp_fst_list,fst_select)

}

head(snp_fst_list)
snp_fst_list<-snp_fst_list %>% arrange(WEIR_AND_COCKERHAM_FST) %>%
  separate_wider_delim(populations, "_", names = c("pair1",NA,"pair2"))
snp_fst_list$pair1<-sub("^pop", "", snp_fst_list$pair1)
snp_fst_list$pair2<-sub("^pop", "", snp_fst_list$pair2)

head(snp_fst_list)

size<-read.table("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/info/body_size.txt",header=T)
head(size)

summary<-size %>% group_by(pop) %>% summarize(mean = mean(body_size)) %>% arrange(mean)
summary

snp_fst_list %>% filter(pair1 == 7 & pair2 == 9)

size<-size %>% mutate(sample = paste0(pop,".",individual))

size %>% group_by(pop) %>%
ggplot(.,aes(x=pop,y=body_size)) + geom_boxplot()

ggplot(size, aes(x=as.factor(pop), y=body_size, fill=as.factor(pop))) +
  geom_boxplot()

head(size)
head(summary)
names(summary)[2]<-"mean_pair2"
joined_fst2<-left_join(joined_fst,summary)
names(joined_fst)[6]<-"mean_pair1"

snp_fst_list$pair1<-as.numeric(snp_fst_list$pair1)
snp_fst_list$pair2<-as.numeric(snp_fst_list$pair2)

head(joined_fst2)

joined_test<-joined_fst2 %>% mutate(ratio = mean_pair1/mean_pair2)
joined_test %>% filter(POS == "3185216") %>%
  ggplot(.,aes(y=WEIR_AND_COCKERHAM_FST,x=ratio)) + geom_point() +geom_smooth(method = "lm")

