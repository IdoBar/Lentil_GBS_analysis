# install needed pacjages/functions

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
gff_file <- "../../../L_culinaris_genome/Lensculinaris_1.2b_genes_fixed_phase.gff3"
gff3 <- read_tsv(gff_file, col_names = c("seqid", "source", "type", "start", "end", "score", "strand",
                                           "phase", "attributes")) 


# gff3 %>% filter(!grepl(paste(discard_pattern, collapse = "|"), attributes))
txdb <- makeTxDbFromGFF(gff3_file, format="gff3", dataSource = "v1.2 UoS",
                        organism = "Lens culinaris",
                        taxonomyId = 3864)
# Match seqlevels between txdb and vcf
# seqlevels(txdb) <- toupper(seqlevels(txdb))
# seqlevels(txdb) <- sub("Arab_Me14_", "", seqlevels(txdb))
# Break the vcf by each isolate

# Find all variants for each isolate and add it up to a list of data frames
stacks_name <- "ref_stacks2"
recent_vcf <- recent_file("./data/intermediate_files/", glue::glue(".*{stacks_name}.+\\.clean.vcf"))

snp_vcf <- VariantAnnotation::readVcf(recent_vcf)
# seqlevels(snp_vcf) <- seqlevels(txdb)
sample_names <- samples(header(snp_vcf))
parents <- sample_names[grepl("ILL.+RF", sample_names)]
# Get a dataframe with all types of variants
var_types <- c("Coding", "Intron", "FiveUTR", "ThreeUTR", "Intergenic", "Promoter")
getAllVariantsDF <- function(vcf, txdb, types){
  # Get all variant types and combine into a large Granges object
  allvars <- map(types, function(typ) locateVariants(vcf, txdb, get(glue::glue("{typ}Variants"))())) %>%
    # lmap(getMethod(c, "GenomicRanges"))
    do.call(getMethod(c, "GenomicRanges"), .)
  
  return(as.data.frame(allvars,row.names = 1:length(allvars)))
}

getAllVariantsGR <- function(vcf, txdb, types){
  # Get all variant types and combine into a large Granges object
  allvars <- map(types, function(typ) locateVariants(vcf, txdb, get(glue::glue("{typ}Variants"))())) %>%
    # lmap(getMethod(c, "GenomicRanges"))
    do.call(getMethod(c, "GenomicRanges"), .)
  
  # return(as.data.frame(allvars,row.names = 1:length(allvars)))
}
vars_by_sample <- DataFrameList(lapply(parents,
                     function(s) getAllVariants(snp_vcf[,s], txdb, var_types))) %>%
  lapply(., S4Vectors::expand, keepEmptyRows = TRUE) %>% DataFrameList()
varsGR <- getAllVariantsGR(snp_vcf[,parents[1]], txdb, var_types)
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
gmap <- readr::read_csv(glue::glue("./data/qtl2_files/Lentil_GBS_{stacks_name}_gmap.csv"))
marker_map <- as.data.frame(rowRanges(snp_vcf)) %>% rownames_to_column("var_id") # %>% as.tibble()
snp_ranges <- rowRanges(snp_vcf)
# Find QTL regions for the trait on interest and manually select range
qtl_trait <- "Leaf_lesion_percent"
qtl_dpi <- "14"
qtl_region_markers <- function(gmap, region_chr, region_start, region_end){
  gmap %>%
    dplyr::filter(chr==region_chr, pos>=region_start, pos<=region_end)
} 
qtl_regions <- LOD_peaks %>% filter(lodcolumn==qtl_trait, dpi==qtl_dpi) %>% 
  .[,c("chr", "ci_lo", "ci_hi")] %>% pmap_dfr(~qtl_region_markers(gmap, ..1, ..2,..3))


qtl_region_snps <- snp_ranges[names(snp_ranges) %in% qtl_regions$marker]
  
# qtl_region <- gmap %>% filter(chr==sel_qtl$chr, pos>=sel_qtl$ci_lo, pos<=sel_qtl$ci_hi)
# boundary_qtl_snps <- c("tag_144852_343285964", "tag_148476_364570997")
# marker_map <-  read_tsv("data/qtl2_files/Lentil_GBS_pmap.csv")
# Load marker map (extracted from gstacks.fa)
# row.names(snp_vcf)
# 
# stacks_dir <- "../Analysis/ref_stacks/stacks2_population_04_12_2017"
# marker_map <- read_tsv(file.path(stacks_dir, "tag_chrom.map"))


qtl_region_tx <- transcriptsByOverlaps(txdb, qtl_region_snps, type = "any", maxgap = 2e4)

#   
#   extractTranscriptSeqs(fa, qtl_region_tx)
#                                   cdsBy(txdb, by="tx", use.names=TRUE))
# 
# 
# qtl_boundaries <-  marker_map %>% 
#   filter(var_id %in% qtl_region_markers$marker) %>% # , seqnames=="LcChr2"
#   arrange(start) %>% droplevels()
# # extract all SNPs between the 2 regions
# surrounding_interval <- 0
# bound_summary <- qtl_boundaries %>% group_by(seqnames) %>% 
#   summarise(qtl_start=min(start), qtl_end=max(end)) %>% mutate(qtl_range=qtl_end-qtl_start)
# for (chr in levels(bound_summary$seqnames)){
#   chr="LcChr2"
#   start_OL <- bound_summary %>% filter(seqnames==chr) %>% .[,"qtl_start"]
# }

qtl_coding <- codingDF %>% filter(varname %in% names(qtl_region_snps)) %>%
  dplyr::select(seqnames, start, end, width, strand, CONSEQUENCE)
qtl_snps <- as.data.frame(varsDF) %>% dplyr::select(-name) %>% #head()
  inner_join(x = data.frame(TXID=as.character(qtl_region_tx$tx_id), tx_name=qtl_region_tx$tx_name), y = .) %>% distinct() %>% 
  left_join(qtl_coding) 
xlsx::write.xlsx(qtl_snps, 
                 filedate(glue::glue("{stacks_name}_QTL_snps"), 
                            ".xlsx", "QTL_results"), row.names = FALSE)


# Extract genomic region under the QTL
# snp_tx <- transcripts(txdb, filter=list(tx_id=unique(qtl_snps$TXID)))
gene_ids <- sub("\\.\\d+$", "", qtl_region_tx$tx_name)
# gene_pattern <- paste(glue::glue('ID={gene_ids};'), collapse = "|")
qtl_genes <- gff3 %>% mutate(attributes=gsub('\"', '', attributes, fixed = TRUE)) %>% 
  separate(attributes, c("ID", "Name", "Description"), sep = ";.+?=") %>%
  filter(type=="gene", Name %in% gene_ids) #%>% 
# Save to excel table
xlsx::write.xlsx(as.data.frame(qtl_genes), 
                 filedate(glue::glue("{stacks_name}_QTL_genes"), 
                          ".xlsx", "QTL_results"), row.names = FALSE)

#### Extract transcripts from genome ####
tx_seqs <- getSeq(fa, qtl_region_tx)
names(tx_seqs) <- qtl_region_tx$tx_name
# cds_seqs <- extractTranscriptSeqs(fa,
#                         cdsBy(txdb, by="tx", use.names=TRUE)) %>%
#   .[snp_tx$tx_name]
# tx_seqs <- cds_seqs[snp_tx$tx_name]


# keep only chromosomes of interest



# Save associated SNPs transcripts to file
# seqinr::write.fasta(as.list(translate(cds_seqs)), names = names(tx_seqs),
#           file.out = filedate(glue::glue("{stacks_name}_QTL_genes_cds"), 
#                                         ".faa", "QTL_results"))
seqinr::write.fasta(as.list(tx_seqs), names = names(tx_seqs),
          file.out = filedate(glue::glue("{stacks_name}_QTL_genes_txs"), 
                                        ".fna", "QTL_results"))
save.image(filedate(glue::glue("Lentil_GBS_{stacks_name}_QTL_SNP_annotation"), ext = ".RData", outdir = "data/intermediate_files"))
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