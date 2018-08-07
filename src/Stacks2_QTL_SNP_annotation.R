# install needed pacjages/functions

devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")

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

# gff3 %>% filter(!grepl(paste(discard_pattern, collapse = "|"), attributes))
txdb <- makeTxDbFromGFF(fixed_gff3_file, format="gff3", dataSource = "v1.2 UoS",
                        organism = "Lens culinaris",
                        taxonomyId = 3864)
# Match seqlevels between txdb and vcf
# seqlevels(txdb) <- toupper(seqlevels(txdb))
# seqlevels(txdb) <- sub("Arab_Me14_", "", seqlevels(txdb))
# Break the vcf by each isolate

# Find all variants for each isolate and add it up to a list of data frames
stacks_name <- "M3m4n3"
recent_vcf <- recent_file("./data/intermediate_files/", glue::glue(".*{stacks_name}.+\\.clean.vcf"))

snp_vcf <- VariantAnnotation::readVcf(recent_vcf)
# seqlevels(snp_vcf) <- seqlevels(txdb)
sample_names <- samples(header(snp_vcf))
parents <- sample_names[grepl("ILL.+RF", sample_names)]
# Get a dataframe with all types of variants
var_types <- c("Coding", "Intron", "FiveUTR", "ThreeUTR", "Intergenic", "Promoter")
getAllVariants <- function(vcf, txdb, types){
  # Get all variant types and combine into a large Granges object
  allvars <- map(types, function(typ) locateVariants(vcf, txdb, get(glue::glue("{typ}Variants"))())) %>%
    # lmap(getMethod(c, "GenomicRanges"))
    do.call(getMethod(c, "GenomicRanges"), .)
  
  return(as.data.frame(allvars,row.names = 1:length(allvars)))
}

vars_by_sample <- DataFrameList(lapply(parents,
                     function(s) getAllVariants(snp_vcf[,s], txdb, var_types))) %>%
  lapply(., S4Vectors::expand, keepEmptyRows = TRUE) %>% DataFrameList()

# test_DF_list <- vars_by_sample %>%
#   lapply(., S4Vectors::expand, keepEmptyRows = TRUE) %>% DataFrameList()
# varsDF <- getAllVariants(snp_vcf[,s], txdb, var_types)
# vardf1 <- S4Vectors::expand(vars_by_sample[[1]], keepEmptyRows = TRUE)
# 
# stack(test_DF_list)
# Combine tables
varsDF <- stack(vars_by_sample)


# genome_seqs <- read.fasta(genome_fa)
coding <- predictCoding(snp_vcf, txdb, fa) # FaFile(genome_fa)

codingDF <- as.data.frame(coding, row.names = 1:length(coding)) %>% mutate(varname=sub("(ctg\\d+):(\\d+)_.+", "\\1.\\2", names(coding)))

#### QTL SNPs ####
# Retrieve SNPs under the QTL (any snps between tag_144852_343285964 to tag_148476_364570997)
LOD_peaks <- readxl::read_excel(recent_file("./QTL_results", 
                                            glue::glue("LOD_peaks_{stacks_name}.+\\.xlsx")))
gmap <- readr::read_csv(glue::glue("./data/qtl2_files/Lentil_GBS_stacks_{stacks_name}_gmap.csv"))
marker_map <- as.data.frame(rowRanges(snp_vcf)) %>% rownames_to_column("var_id") %>% as.tibble()
sel_qtl <- LOD_peaks %>% filter(lod==max(lod))
qtl_region <- gmap %>% filter(chr==sel_qtl$chr, pos>=sel_qtl$ci_lo, pos<=sel_qtl$ci_hi)
# boundary_qtl_snps <- c("tag_144852_343285964", "tag_148476_364570997")
# marker_map <-  read_tsv("data/qtl2_files/Lentil_GBS_pmap.csv")
# Load marker map (extracted from gstacks.fa)
# row.names(snp_vcf)
# 
# stacks_dir <- "../Analysis/ref_stacks/stacks2_population_04_12_2017"
# marker_map <- read_tsv(file.path(stacks_dir, "tag_chrom.map"))
qtl_boundaries <-  marker_map %>% filter(var_id %in% qtl_region$marker) %>% # , seqnames=="LcChr2"
  arrange(start) %>% droplevels()
# extract all SNPs between the 2 regions
surrounding_interval <- 2e6
bound_summary <- qtl_boundaries %>% group_by(seqnames) %>% 
  summarise(qtl_start=min(start), qtl_end=max(end)) %>% mutate(qtl_range=qtl_end-qtl_start)
for (chr in levels(bound_summary$seqnames)){
  chr="LcChr2"
  start_OL <- bound_summary %>% filter(seqnames==chr) %>% .[,"qtl_start"]
}


qtl_snps <- as.data.frame(varsDF) %>% dplyr::select(-name) %>%
  # mutate_if(is.list, possibly(`[[`, NA)) %>%
  filter(seqnames %in% qtl_boundaries$seqnames,
                              start >= qtl_boundaries$start[1], 
                              start <= qtl_boundaries$start[nrow(qtl_boundaries)]) %>% distinct()
# Extract genomic region under the QTL
snp_tx <- transcripts(txdb, filter=list(tx_id=unique(qtl_snps$TXID)))
gene_ids <- sub("\\.\\d+$", "", snp_tx$tx_name)
# gene_pattern <- paste(glue::glue('ID={gene_ids};'), collapse = "|")
#### Extract transcripts from genome ####
cds_seqs <- extractTranscriptSeqs(fa,
                                  cdsBy(txdb, by="tx", use.names=TRUE))
tx_seqs <- cds_seqs[snp_tx$tx_name]
snp_genes <- gff3 %>% mutate(attributes=gsub('\"', '', attributes, fixed = TRUE)) %>% 
  separate(attributes, c("ID", "Name", "Description"), sep = ";.+?=") %>%
  filter(type=="gene", Name %in% gene_ids) #%>% 
  
  # mutate(ID=sub("ID=(.+?);.+", "\\1", attributes)) %>% 
  

head(codingDF)


save.image(filedate(glue::glue("Lentil_GBS_{stacks_name}_QTL_SNP_annotation"), 
                    ext = ".RData", outdir = "data"))



# keep only chromosomes of interest



# Save associated SNPs transcripts to file
seqinr::write.fasta(as.list(translate(tx_seqs)), names = names(tx_seqs),
                    file.out = filedate(glue::glue("{stacks_name}_QTL_genes_cds"), 
                                        ".fasta", "QTL_results"))

#### Prepare summary table ####
# Load BLAST results
# blast_fields <- unlist(strsplit("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore subj_desc subj_sci_name subj_kingdom", split=" "))  # pcov
# annot_table <- read_tsv("data/Assoc_SNPs/Assoc_SNP_genes_cds.nr_blastp.outfmt6",
#                 col_names = blast_fields)
loci_annotation <- readxl::read_excel('data/Assoc_SNPs/Loci_annotation.xlsx', "Tx_annotation") %>%
  group_by(qseqid) %>% summarise(gene_description=paste(subj_desc, collapse="|")) %>%
  left_join(as.data.frame(assoc_snp_tx), ., by=c("tx_name"="qseqid")) %>%
  mutate(gene_id=as.character(gene_id))
# columns(txdb)
assoc_table <- left_join(assoc_snp_annot, loci_annotation[c("gene_id", "gene_description")]) #  %>% 
# dplyr::select(-name) %>% as_tibble() %>% distinct()       #  by=c("GENEID"="gene_id"))
write_tsv(assoc_table, filedate("Assoc_SNP_table_annotated", ".txt", "data/Assoc_SNPs"))