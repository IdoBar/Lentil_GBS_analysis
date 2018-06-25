devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# Install and load required packages from GitHub
# .fixdevtools()
# remove.packages("onemap")
options(unzip = "/usr/bin/unzip")
Sys.setenv(TAR = "/bin/tar")
# pacman::p_load_gh(c("IdoBar/onemap@patch-3")) # "thierrygosselin/radiator", "augusto-garcia/onemap"
# vignette("Inbred_Based_Populations") # OneMap tutorial
pacman::p_load_gh("augusto-garcia/onemap") # devtools::install_github
# pacman::p_load(Rhtslib)
# Install and load needed packages
package_list <- c("tidyverse", "qtl", "ASMap", "RColorBrewer", "Rsamtools", "doFuture", 
                  "foreach", "plotly") # "radiator", 
pacman::p_load(char=package_list)

# Onemap requires R>=3.4 and Rhtslib and zlibbioc - install with conda
# git_packs <- c("augusto-garcia/onemap")
# install.deps(git_packs, repo = "git")
# devtools::install_github(git_packs)

options(stringsAsFactors=FALSE)

# Load vcf file
stacks_dir <- "../Stacks_vcf_files"
vcf_file <- "ref_stacks2_populations.snps.vcf"  #
vcf_data <- read_tsv(file.path(stacks_dir, vcf_file), comment = "##") #%>% 
  # mutate(`#CHROM`=as.numeric(unlist(strsplit(`#CHROM`, split="_"))[2]))  # n_max=100
sample_cols <- colnames(vcf_data)[10:ncol(vcf_data)]
# Relace missing genotypes and heterozygotes with ./.
vcf_data <- vcf_data %>% mutate_at(vars(one_of(sample_cols)), 
       .funs = funs(if_else(grepl("^1/0", .) | grepl("^0/1", .) | grepl("^\\.$", .), "./.", .)))

# Check if we have any heterozygote left
# sum(apply(vcf_data[sample_cols], 1, function(r) sum(grepl("^1/0", r) || grepl("^0/1", r))))
exclude_markers <- FALSE
vcf_filt_data <- vcf_data
incl_samples <- sample_cols
counter=1
miss_rates <- seq(0.2, 0.1, -0.05)
for (m in miss_rates){
  
  # sample_col_filt <- colnames(vcf_filt_data)[10:ncol(vcf_filt_data)]
  # m <- miss_rates[3]
  
  excl_samples <- incl_samples
  
  while (length(excl_samples)>0) {
    # vcf_filt_data <- vcf_filt_data %>% filter(!exclude_markers)
    exclude_markers <- apply(vcf_filt_data[incl_samples], 1, function(marker) {
      sum(grepl("^\\./\\.", marker) | grepl("^\\.$", marker)) > m*length(incl_samples)
    } )
    table(exclude_markers)
    
    vcf_filt_data <- vcf_filt_data %>% filter(!exclude_markers)
    # Find samples with too many missing markers
    geno_rate <- 0.8
    excl_samples <- incl_samples[apply(vcf_filt_data[incl_samples], 2, 
                                       function(geno) sum(grepl("\\./\\.", geno))>(1-geno_rate)*nrow(vcf_filt_data))]
    
    
    # Now repeat marker exclusion without the problematic genotypes
    incl_samples <- incl_samples[!incl_samples %in% excl_samples]
    LogMsg(sprintf("Genotype and Marker filtering, round %d, missing rate %.2f", counter, m))
    LogMsg(sprintf("Filtered %d genotypes and %d markers", length(excl_samples), sum(exclude_markers)))
    counter <- counter + 1
  }
  
}
parents <- incl_samples[grepl("_RF", incl_samples)]
poly_markers <- apply(vcf_filt_data[parents], 1, function(g) {
  # g <- vcf_data[1, parents]
  genotypes <- sub(":.+", "", g)
  # log_vec <- ifelse(biallelic, length(prop_names[!grepl("\\.", prop_names)])==2,
  # length(prop_names[!grepl("\\.", prop_names)])>1)
  p_table <- prop.table(table(genotypes))
  g_table <- p_table[!grepl("\\./\\.", names(p_table))]
  if (length(g_table)==2) {
    res <- TRUE # g_table[1]>0.05 && g_table[2]>0.05
  }  else res <- FALSE
  return(res)
})

table(poly_markers)

# vcf_filt_data <- vcf_filt_data %>% filter(poly_markers)
exclude_genotypes <- sample_cols[!sample_cols %in% incl_samples]
# Remove missing genotypes and markers, update marker information (CHROM, POS, LOC)
clean_vcf <- vcf_filt_data %>% filter(poly_markers) %>% dplyr::select(-one_of(exclude_genotypes)) #%>% 
  # mutate(ID=paste(ID, POS, sep= "_"))

# Load marker map (extracted from gstacks.fa)
# marker_map <- read_tsv(file.path(stacks_dir,"tag_chrom.map")) %>% 
#   dplyr::rename(CHROM_POS=POS) 
  # left_join(marker_map, c("#CHROM" = "LOC")) %>% 
  # mutate(ID=paste0("tag",`#CHROM`, POS, sep="_"), `#CHROM` = CHROM, 
  #        POS=ifelse(STRAND=="+", CHROM_POS+POS, CHROM_POS-POS)) %>% 
  # dplyr::select(-one_of(colnames(marker_map)))

# Write back vcf
vcf_basename <- tools::file_path_sans_ext(vcf_file)
clean_vcf_file <- filedate(sprintf("%s.miss%d",vcf_basename, m*100), ".clean.vcf","data", 
                           dateformat = FALSE)
                            # sub(".clean.sorted.vcf", "_miss40.clean.vcf",vcf_file, fixed = TRUE))

# unlink(clean_vcf_file)
# Copy the header into a new file
system2("grep", args=c("'^##'", file.path(stacks_dir, vcf_file)), stdout = clean_vcf_file)
if (file.exists(clean_vcf_file)) write_tsv(clean_vcf, clean_vcf_file, append = TRUE, col_names = TRUE)


# Import to OneMap 
# Try to read straight from vcf
clean_vcf_file <- "data/ref_stacks2_populations.snps.miss10.clean.vcf"
clean_vcf_basename <- sub(".clean.vcf", "", basename(clean_vcf_file))
onemap_file <- onemap_read_vcfR(vcfR::read.vcfR(clean_vcf_file), cross = "ri self",
                                parent1 = "ILL6002_RF", 
                                parent2 = "ILL7537_RF")
# rm(clean_vcf, clean_vcfR, vcf_data, vcf_geno)
# Try to convert into MapMaker format 
# vcf2raw(input = paste0(clean_vcf_file, ".bgz"), output = paste0(clean_vcf_file, ".raw"), 
#         cross = "ri self",
#         parent1 = "ILL6002_RF", 
#         parent2 = "ILL7537_RF")
# 
# onemap_file <- read_onemap(dir=dirname(clean_vcf_file), 
#                            paste0(basename(clean_vcf_file), ".raw"))
# # Try to load the file produced by Vcf2MapMaker.pl
# onemap_file <- read_mapmaker(dir=dirname(vcf_file), "gstacks_no_RL_NS_GQ_DP_Samples_filtered.mapmaker.raw" )#              sub(".vcf", ".mapmaker.raw", basename(vcf_file)))
# "gstacks_case4_Hom_Pol.imputed.clean.mapmaker.raw"
# class(onemap_file)

# plot the data
# plot.onemap(onemap_file)
# plot_by_segreg_type(onemap_file)
# Check for marker segregation
ril_test <- test_segregation(onemap_file) 
Bonferroni_alpha(ril_test)
plot(ril_test)
good_markers <- select_segreg(ril_test, distorted = FALSE, numbers = TRUE)

(LOD_sug <- suggest_lod(onemap_file)) # estimate LOD threshold
twopts_ril <- rf_2pts(onemap_file, LOD=floor(LOD_sug), max.rf = 0.25)
# Assigning markers to linkage groups
mark_all_ril <- make_seq(twopts_ril, good_markers)
# LGs_ril <- group(mark_all_ril, max.rf = 0.25)
# LGs_ril <- group(mark_all_ril)

# test a range of LOD and rf values (in parallel)
# Find out how many cores are available (if you don't already know)
# detectCores()
# Create cluster with desired number of cores
# cl <- makeForkCluster(detectCores()/2)
# Register cluster
# registerDoParallel(cl)
# registerDoSEQ()
# Create cluster with desired number of cores
registerDoFuture()
plan(multisession, workers=8)

LOD_range <- 4:10
rf_range <- seq(0.1,0.5, 0.05)
LogMsg(sprintf("Starting testing LOD and max.rf parameters.", paste(LOD_range, collapse = ","),
               paste(rf_range, collapse = ",")))
LGs_df <- foreach(l = LOD_range, .combine = bind_rows) %:% 
  foreach(rf = rf_range, .combine = bind_rows) %dopar% {
    # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
    # l=3, rf=0.2
    dummy <- length(onemap_file)
    LG <- onemap::group(mark_all_ril, max.rf = rf, LOD = l)
    data.frame(LOD=l, max.rf=rf, as.data.frame(table(LG$groups)))
    
  }
LogMsg(sprintf("Finished testing LOD (%s) and max.rf (%s) parameters.", paste(LOD_range, collapse = ","),
               paste(rf_range, collapse = ",")))
  # LGs_df <- bind_rows(LGs_df, df)
# stopCluster(cl)


colnames(LGs_df) <- c("LOD", "max.rf", "LG", "Freq")
write_tsv(LGs_df, filedate(sprintf("%s_LOD_rf_params",clean_vcf_basename), ".tsv", "data/itermediate_files"))
LGs_df <- read_tsv("data/intermediate_files/ref_stacks2_populations.snps.miss10_LOD_rf_params_13_06_2018.tsv")
LGs_sum <- LGs_df  %>% group_by(LOD, max.rf) %>% summarise(n=n()) %>% mutate(LOD_cat=as.character(LOD))

LGs_df %>% filter(LOD>7)
p <- ggplot(LGs_sum, aes(x=max.rf, y=n, color=LOD_cat)) + geom_line(lwd=1) + scale_color_brewer(palette = "Paired")
ggplotly(p)

l=4 # 6
rf=0.35  #0.25

set_map_fun(type = "kosambi") # haldane
LG <- group(mark_all_ril, max.rf = rf, LOD = l)

# Check markers distribution in Linkage groups
table(LG$groups)
# Order the markers in first LG
# LG1_f2 <- make_seq(LG, 1)
# # Construct the linkage map, by automatic usage of the try algorithm:
# LG1_f2_ord <- order_seq(input.seq = LG1_f2, n.init = 5,
#                         subset.search = "twopt",
#                         twopt.alg = "rcd", THRES = 3,
#                         touchdown = TRUE)
# Get the order with all markers:
# (LG1_f2_final <- make_seq(LG1_f2_ord, "force"))
# Make a list to hold all the linkage groups
LG_list <- list()
# g=1
# 
# LG_list[[paste0("LG", g)]] <- LG1_f2_final
# ripple_seq(LG1_f2_final, ws = 5, LOD = 3)
# Do the same for the rest of the linkage groups
for (g in 1:LG$n.groups){
  # Order the markers in the current LG
  LGg_f2 <- make_seq(LG, g)
  # Construct the linkage map, by automatic usage of the try algorithm:
  LogMsg(sprintf("Determining order of LG %d, please wait", g))
  LGg_f2_ord <- order_seq(input.seq = LGg_f2, n.init = 5,
                          subset.search = "twopt",
                          twopt.alg = "rcd", THRES = 3,
                          touchdown = TRUE)
  LogMsg(sprintf("Finished ordering LG %d!", g))
  # Get the order with all markers:
  LG_list[[paste0("LG", g)]] <- make_seq(LGg_f2_ord, "force")
}

save.image(filedate("Lentil_GBS_linkage_map", ".RData", "data"))

load('data/Lentil_GBS_linkage_map_02_03_2018.RData')

# View linkage map
draw_map(LG_list, names = TRUE, grid = FALSE, cex.mrk = 0.7)
# Number of markers in each Chromosome:
sapply(LG_list, function(lg) length(lg$seq.num))

# match chromosome positions of markers:
markers <- colnames(onemap_file$geno)
temp_map_file <- tempfile()
write_map(LG_list,temp_map_file)
temp_map <- read_delim(temp_map_file, delim = " ",
           col_names = c("LG", "marker", "pos")) %>%
  mutate(Phys_Chrom=onemap_file$CHROM[match(marker, markers)], 
         Phys_Pos=onemap_file$POS[match(marker, markers)], marker_num=match(marker, markers)) %>% 
  write_csv(filedate("Temporary_linkage_map_order", ".csv", "data"))


# Manually adjust markers at edges of LGs:
# View the LOD scores and RF of LG4
rf_graph_table(LG_list$LG4)
LG_list$LG4 <- drop_marker(LG_list$LG4, c(960:967))
LG4_ord <- order_seq(input.seq = LG_list$LG4, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)
LogMsg(sprintf("Finished ordering LG %d!", g))
# Get the order with all markers:
LG4_seq <- make_seq(LG4_ord, "force")
# Checking for better orders: (seems like the mapping order is pretty good)
# ripple_seq(LG4_seq, ws = 5)
# # try to add them again
# LG4_add <- try_seq(LG4_seq, 967)
# Insert back into LG_list
LG_list$LG4 <- LG4_seq

# View the LOD scores and RF of LG2
rf_graph_table(LG_list$LG2)
LG_2_trim <- drop_marker(LG_list$LG2, c(438:439, 531:534))
LG2_ord <- order_seq(input.seq = LG_2_trim, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "rcd", THRES = 3,
                     touchdown = TRUE)
LogMsg(sprintf("Finished ordering LG %d!", g))
# Get the order with all markers:
LG2_seq <- make_seq(LG2_ord, "force")
# Insert back into LG_list
LG_list$LG2 <- LG2_seq
# View linkage map
draw_map(LG_list, names = FALSE, grid = FALSE, cex.mrk = 0.7)

#### Save files for QTL analysis ####
# Export linkage map
write_map(LG_list,filedate("Lentil_GBS_linkage_map_order", ".map", "data", dateformat = FALSE))






#### Chromosome-based mapping ####
Chroms <- unique(twopts_ril$CHROM)
CHR_mks <- group_seq(input.2pts = twopts_ril, # seqs = "CHROM", unlink.mks = mark_all_ril,
                     rm.repeated = FALSE, LOD = l, max.rf = rf)




# Make a list to hold all the linkage groups
Chr_list <- list()
Chroms <- unique(twopts_ril$CHROM)
# ripple_seq(LG1_f2_final, ws = 5, LOD = 3)
# Do the same for the rest of the linkage groups
for (chr in Chroms){
  # Order the markers in the current LG
  # chr=Chroms[1]
  # Chr_f2 <- make_seq(CHR_mks$sequences[[paste0("CHR",chr)]])
  # Construct the linkage map, by automatic usage of the try algorithm:
  LogMsg(sprintf("Determining order of Chrom %s, please wait", chr))
  Chr_ord <- order_seq(input.seq = CHR_mks$sequences[[paste0("CHR",chr)]], n.init = 5,
                          subset.search = "twopt",
                          twopt.alg = "rcd", THRES = 3,
                          touchdown = TRUE)
  LogMsg(sprintf("Finished ordering Chrom %s!", chr))
  # Get the order with all markers:
  Chr_list[[chr]] <- make_seq(Chr_ord, "force")
}

# View the LOD scores and RF
rf_graph_table(Chr_list$LcChr2, inter = FALSE)

draw_map(Chr_list, names = FALSE, grid = FALSE, cex.mrk = 1)

# Number of markers in each Chromosome:
sapply(LG_list, function(lg) length(lg$seq.num))

#### ASMap ####
# cleaned file to read
clean_vcf_file <- "data/ref_stacks2_populations.snps.miss15.clean.vcf"
# Manually convert to MapMaker format
clean_vcf <- read_tsv(clean_vcf_file, comment = "##")
# Subset just the genotypes (also remove unwanted samples, RL and genotype calling meta-information)
vcf_geno <-  clean_vcf[10:ncol(clean_vcf)] %>% 
  mutate_all(.funs = funs(sub(":.+", "", .))) # clean_vcf[incl_samples]

# Create the dictionary to replace the values
dict_keys <- c("0/0", "0/1", "1/0", "1/1", "./.")
# For ASMap
mst_dict <- c("A", "X", "X", "B", "U")
# For MapMaker
MapMaker_dict <- c("A", "H", "H", "B", "-")

# geno_converted <- dict_values[geno]

# Apply the conversion to the entire table
# create a custom function to perform the conversion on a vector of input strings
convert_genotypes <- function(genos, diction, keys){
  names(diction) <- keys
  # Convert the 1st object in each list (the genotype) based on our dictionary
  geno_converted <- sapply(genos, function(g) {
    diction[sub(":.+", "", g)]
  })
  return(geno_converted) # return a vector of converted genotypes
}

# Import to ASMap
raw_geno <- as.data.frame(apply(vcf_geno, 2, function(col) convert_genotypes(col, mst_dict,
                                                  dict_keys)), stringsAsFactors=FALSE)
# row.names(raw_geno) <- paste(clean_vcf$`#CHROM`, clean_vcf$POS, sep="_") 
row.names(raw_geno) <- clean_vcf$ID
raw_geno[1:10,1:10]


testd <- mstmap(raw_geno, pop.type = "RIL5", miss.thresh=0.4, p.value = 5e-3) # p.value = 5e-2

gc <- genClones(testd, tol = 0.95)
if ("cgd" %in% names(gc)) {
  mapBC2 <- fixClones(testd, gc$cgd, consensus = TRUE)
} else {mapBC2 <- testd}
profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Remove highly distorted markers
mm <- statMark(mapBC2, stat.type = "marker")$marker$AB
dm <- markernames(mapBC2)[(mm > 0.9) | (mm < 0.1)]
mapBC3 <- drop.markers(mapBC2, dm)
mapBC4 <- jittermap(mapBC3)
profileMark(mapBC4, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
mapBC5 <- mstmap(mapBC4, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = 1e-3)



