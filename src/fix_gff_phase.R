devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")

if(getRversion() >= "3.5.0") {
  BiocManager::install(c("GenomicFeatures", "VariantAnnotation", "GenomicRanges", "Rsamtools", "plyranges"))
  pacman::p_load(plyranges)
}


pacman::p_load(char = c( "seqinr", "tidyverse", "GenomicFeatures", "VariantAnnotation", "GenomicRanges", "Rsamtools")) #  "plyranges", "BiocManager" requires R>=3.5
# package.version("GenomicFeatures")
# bioc_packages <- c("GenomicFeatures", "VariantAnnotation", "GenomicRanges", "Rsamtools")
# # pacman::p_load(char = bioc_packages)
# install.deps(bioc_packages, repo = "bioc")

#### Load genome files  ####
# Load genome file
genome_fa <- "../../../L_culinaris_genome/Lensculinaris_genome_v1.2.fasta"
if (!file.exists(paste0(genome_fa, ".fai"))) indexFa(genome_fa)
fa <- open(FaFile(genome_fa))
# load gene models file
gff_file <- "../../../L_culinaris_genome/Lensculinaris_1.2b_genes.gff3"
gff3 <- read_tsv(gff_file, col_names = c("seqid", "source", "type", "start", "end", "score", "strand",
                                         "phase", "attributes")) 
#### Fix missing phasing information ####
missing_phase <- gff3 %>% filter(type=="CDS", phase==".", (end-start)>=18)  %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# mutate(orig_id=seqid, seqid=gsub("ID=(.+);.+", "\\1", attributes)) %>% makeGRangesFromDataFrame()
# system2("gunzip", args = c("-cd", genome_fa), stdout = "../../../L_culinaris_genome/Lensculinaris_genome_v1.2.fasta")
missing_cds <- getSeq(fa, missing_phase)
names(missing_cds) <- missing_phase$attributes
# remove short ones
# missing_cds
# test_seq <- missing_cds[length(missing_cds)]
all_cds_df <- lmap(missing_cds,  ~map(1:3, function(pos) {
  tryCatch({
    tibble(AA_seq=toString(Biostrings::translate(subseq(., start=pos))), 
           phase=pos-1, attributes=names(.))},
    warning = function(c) {warning(c); NULL},
    error = function(c) {warning(c); NULL}
  )
}))  %>% map_df(bind_rows) # %>% filter(!grepl("\\*\\w+", AA_seq))

good_cds <- all_cds_df %>% filter(!grepl("\\*\\w", AA_seq))

# discard_cds <- missing_phase$attributes[!missing_phase$attributes %in% good_cds$attributes]
# discard_pattern <- sub("cds(\\d+;)", "+[a-z]\\1", discard_cds)  # remove just cds and exons
# discard_pattern <- sub("\\.cds(\\d+;).+", "", discard_cds) # remove entire ID
# fix phase information
gff3$phase[gff3$attributes %in% good_cds$attributes] <- good_cds$phase
# Temporary fix (change all unfixable phases to 0)
gff3$phase[gff3$type=="CDS" & gff3$phase=="."] <- 0
fixed_gff3_file <- sub("(\\.gff[3]*)$", "_fixed_phase\\1", gff_file)
gff3 <- gff3 %>% filter(!(type=="CDS" & phase==".")) %>% write_tsv(fixed_gff3_file)