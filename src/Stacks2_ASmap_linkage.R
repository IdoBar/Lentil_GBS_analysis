devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R") # needs Rcurl installed
# Install and load required packages from GitHub
# .fixdevtools()
# remove.packages("onemap")
# devtools::install_github("augusto-garcia/onemap")
# pacman::p_load_gh(c("augusto-garcia/onemap")) # "thierrygosselin/radiator", "augusto-garcia/onemap", "IdoBar/onemap@patch-3"
# vignette("Inbred_Based_Populations") # OneMap tutorial
# devtools::install_github("augusto-garcia/onemap")
# Install and load needed packages
package_list <- c("tidyverse", "qtl", "RColorBrewer", "Rsamtools", "doFuture", "foreach", "ASMap") # "radiator", 
pacman::p_load(char=package_list)

# Onemap requires R>=3.4 and Rhtslib and zlibbioc
# git_packs <- c("augusto-garcia/onemap")
# install.deps(git_packs, repo = "git")
# devtools::install_github(git_packs)

options(stringsAsFactors=FALSE)

# Load vcf file and load associated .RData
stacks_name <- "ref_stacks2"
vcf_file <- recent_file("../Stacks_vcf_files/filtered_VCFs", 
                        glue::glue("{stacks_name}_populations.+.vcf"))  # ref_stacks2_populations.snps.vcf
group_markers <- "LG"
vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))
  
source("./Stacks2_linkage_mapping/src/vcf_filtration.R")
vcf_filtration(vcf_file, miss_rates = seq(0.2, 0.15, -0.05), geno_rate = 0.8)

#### ASMap ####
# cleaned file to read
# clean_vcf_file <- "./data/intermediate_files/ref_stacks2_populations.snps.miss10.clean.vcf"
clean_vcf_basename <- sub(".clean.vcf", "", basename(clean_vcf_file))
# stacks_name <- sub("_populations.snps", "", sub("\\.miss.+", "", clean_vcf_basename))
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

range_length <- 25
# p_value_range <- 5e-10*(5^(1:range_length))
p_value_range <- 1e-12*(2^(0:range_length-1))
# p_value_range <- 1e-8*(2^(0:(range_length-1)))

# Create cluster with desired number of cores
registerDoFuture()
plan(multiprocess, workers=3)

# asmap_thres_df <- foreach(p = p_value_range, .combine = bind_rows) %dopar% {
#     # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
#     # l=3, rf=0.2
#     # dummy <- length(onemap_file)
#     testd <- R.utils::captureOutput(mstmap(raw_geno, pop.type = "RIL5", miss.thresh=0.1, p.value = p))
#     # p.value = 5e-2
#     
#     data.frame(p.value = p, linkage_groups=as.numeric(stringr::str_extract(testd[1], "\\d+$")))
#     
#   }

asmap_thres_list <- foreach(p = p_value_range, 
          .final = function(x) setNames(x, paste0("asmap_p", seq_along(p_value_range)))) %dopar% {
  # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
  # l=3, rf=0.2
  # dummy <- length(onemap_file)
  mstmap(raw_geno, pop.type = "RIL5", miss.thresh=0.1, p.value = p, bychr = FALSE) 
  # p.value = 5e-2 , objective.fun = "ML"
  
}

# Check number of markers per groups
sapply(asmap_thres_list, function(l) qtl::nmar(l))
prettyNum(p_value_range)
# process_asmap <- function(asmap){
#   # asmap <- asmap_thres_list[1]
#   tibble(map_name=names(asmap), linkage_group=names(qtl::nmar(asmap[[1]])), 
#              marker_num=qtl::nmar(asmap[[1]]))
# }

asmap_thres_df <- imap_dfr(asmap_thres_list, ~tibble(map_name = .y, linkage_group = names(qtl::nmar(.x)), 
                                                     marker_num = qtl::nmar(.x)))

asmap_thres_df %>% filter(marker_num>150) %>% count(map_name)

# asmap_thres_list %>% map_df( ~ as.list(qtl::nmar(.x[[1]]))) 
# 
# asmap_thres_df[2,2]
# 
# ~ tibble(p_value = names(.x), linkage_group=names(qtl::nmar(asmap_thres_list[[p]])), marker_num=qtl::nmar(asmap_thres_list[[p]]))) %>% purrr::map_df(bind_rows)
# paste(qtl::nmar(asmap_thres_list[[p]]), collapse = ";")
# # Check length of groups
# sapply(asmap_thres_list, function(l) qtl::chrlen(l))


# selected p7 (from range)
selected_p <- 20
asmap <- asmap_thres_list[[selected_p]]
LogMsg(glue::glue("Selected p.value is: {p_value_range[selected_p]}"))
LogMsg(glue::glue("Number of markers in each linkage group: {paste(qtl::nmar(asmap), collapse=', ')}"))
LogMsg(glue::glue("Length of each linkage group: {paste(qtl::chrlen(asmap), collapse=', ')}")) 


# use function quickEst() to re-estimate the genetic distances 
qtl::plotMissing(asmap)
sg <- statGen(asmap, bychr = FALSE, stat.type = "miss")
mapBC1 <- subset(asmap, ind = sg$miss < sum(qtl::nmar(asmap))*0.15)
gc <- genClones(mapBC1, tol = 0.975)
if ("cgd" %in% names(gc)) {
  mapBC2 <- fixClones(mapBC1, gc$cgd, consensus = FALSE)
} else {mapBC2 <- mapBC1}
mapBC3 <- jittermap(mapBC2)
# Check 
profileMark(mapBC3, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Remove highly distorted markers
mm <- statMark(mapBC3, stat.type = "marker")$marker$AA
LogMsg(glue::glue("Number of distorted markers:: {sum(mm > 0.7 | mm < 0.3)}"))
dm <- markernames(mapBC3)[(mm > 0.7) | (mm < 0.3)]
mapBC4 <- drop.markers(mapBC3, dm)

profileMark(mapBC4, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Pull markers that are co-located
# mapDHs <- pullCross(mapBC4, type = "co.located")
mapDHs <- mapBC4
# Pull markers that have significant segregation distortion with a p-value less than 0.02 and have a missing
# value proportion greater than 0.03.
mapDHs <- pullCross(mapDHs, type = "seg.distortion", pars =
                      list(seg.ratio = "40:0:60"))
mapDHs <- pullCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.075))
names(mapDHs)
ncol(mapDHs$seg.distortion$data)
ncol(mapDHs$missing$data)
# profileMark(mapDHs, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
nmar(mapDHs)
# Visualise data
heatMap(mapDHs, lmax = 20)
# Profile the statistics for the genotypes across
pg <- profileGen(mapDHs, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id =
             "Genotype", xo.lambda = 15, layout = c(1, 3), lty = 2, cex=0.7)
miss_ind <- pg$stat$miss[pg$stat$miss>400]

mapBC5 <- subsetCross(mapDHs, ind = !pg$stat$miss>400)
mapRIL5 <- mstmap(mapBC5, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = p_value_range[selected_p])


# calculate genetic distance
# map2 <- quickEst(mapDH, map.function = "kosambi")

# Visualise data
lod_range <- seq(10,100, by = 5)

foreach(l=lod_range) %dopar% {
  pdf(filedate(sprintf("%s_asmap_lod%d_test", stacks_name, l), ".pdf", 
               glue::glue("./plots/{stacks_name}/param_estimation")), width = 8, height = 8)
  heatMap(mapRIL5, lmax = l)
  dev.off()
  return(NULL)
  
}
profileMark(mapRIL5, stat.type = c("seg.dist", "prop", "dxo", "recomb"),
            layout = c(1, 5), type = "l")
heatMap(mapRIL5, lmax = 30, chr = "L.13")
nmar(mapRIL5)
chrlen(mapRIL5)
mapRIL5$geno$L.4
# Split linkage group and merge with another (based on heatmap)
mapRIL5_split <- breakCross(mapRIL5, split = list(`L.4` = "305093:50:+"), suffix = "alpha", sep = "")
mapRIL5_merged <- mergeCross(mapRIL5_split, merge = list(`L.4` = c("L.4", "L.13"), 
                             `L.1` = c("L.3","L.2", "L.1")))
# mapRIL5_merged <- mstmap(mapRIL5_merged, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = p_value_range[selected_p])
# Remove markers from L.10,6,7,15
heatMap(mapRIL5_merged, lmax = 30)
nmar(mapRIL5_merged)
chrlen(mapRIL5_merged)
heatMap(mapRIL5_merged, lmax = 30, chr = "L.4")
# Push back markers with higher missing rates
mapRIL5_added <- pushCross(mapRIL5_merged, type = "missing", pars = list(miss.thresh =
                                                                           0.1, max.rf = 0.3))
# visualise
heatMap(mapRIL5_added, lmax = 30)
# reconstruct the map (but don't reassess linkage groups - p.value=2)
mapRIL5b <- mstmap(mapRIL5_added, bychr = TRUE, trace = TRUE, dist.fun = "kosambi",
                   p.value = 2)
mapRIL5b$geno$L.4$map[350:length(mapRIL5b$geno$L.4$map)]
mapRIL5b$geno$L.4$map[500:520]

mapRIL5b_split <- breakCross(mapRIL5b, split = list(`L.4` = c("723124:89:-", "186649:35:-")), suffix = "alpha", sep = "")
mapRIL5b_merged <- mergeCross(mapRIL5b_split, merge = list(`L.4` = c("L.4A", "L.4C"), 
                                                   `L.5` = c("L.4B","L.5")))
mapRIL5b <- mstmap(mapRIL5b_merged, bychr = TRUE, trace = TRUE, dist.fun = "kosambi",
                   p.value = 2)
# Check for problematic genotypes
pg1 <- profileGen(mapRIL5b, bychr = FALSE, stat.type = c("xo", "dxo","miss"), 
                  id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

# remove problematic genotypes
mapRIL5c <- mstmap(subsetCross(mapRIL5b, ind = !pg1$xo.lambda), bychr = TRUE, trace = TRUE, 
                   dist.fun = "kosambi", p.value = 2)
profileMark(mapRIL5c, stat.type = c("seg.dist", "prop", "dxo", "recomb"),
            layout = c(1, 5), type = "l")
# Push back distorted markers
mapRIL5d <- pushCross(mapRIL5c, type = "seg.distortion", pars =
                       list(seg.ratio = "70:0:30"))
mapRIL5e <- mstmap(mapRIL5d, bychr = TRUE, trace = TRUE, dist.fun =
                    "kosambi", p.value = 2)
mapRIL5e_merged <- mergeCross(mapRIL5d, merge = list(`L.1` = c("L.1", "UL2")))
mapRIL5e <- mstmap(mapRIL5e_merged, bychr = TRUE, trace = TRUE, dist.fun =
                     "kosambi", p.value = 2)
# Check difference in linkage group lengths
round(chrlen(mapRIL5e) - chrlen(mapRIL5b), 5)
# Check difference in number of markers 
nmar(mapRIL5e) - nmar(mapRIL5b)
# visualise
heatMap(mapRIL5e, lmax = 30)

#### Compare linkage map to physical map ####
unlinked_markers <- sapply(mapRIL5e$geno, function(lg) {if(length(lg$map>10)) return(names(lg$map))})
marker_map_df <- imap_dfr(mapRIL5e$geno, ~tibble(LG_name = .y, marker_name = markernames(mapRIL5e, .y))) %>%
  inner_join(clean_vcf[c(1:3)], by = c("marker_name" = "ID")) %>% 
  mutate(Chr=sub("_.+", "", `#CHROM`))

marker_map_sum <- marker_map_df %>% group_by(Chr) %>% summarise(num=n())
# Manually merge all linkage groups belonging to the same chromosome
mapRIL5_chrom <- mergeCross(mapRIL5e, merge = list("LcChr1" = c("L.1", "L.8", "L.13"), 
                           "LcChr2" = c("L.10", "L.5", "L.14", "L.15", "L.16"), "LcChr3"=c("L.11", "L.18"),
                           "LcChr4" = c("L.7", "L.19" ), "LcChr5"=c("L.4", "L.21"), 
                           "LcChr6" = c("L.2", "L.12", "L.22"), "LcChr7" = c("L.3", "L.23")))
# visualise
heatMap(mapRIL5_chrom_final, lmax = 30)
chrom_map_df <- imap_dfr(mapRIL5_chrom$geno, ~tibble(LG_name = .y, 
                                     marker_name = markernames(mapRIL5_chrom, .y))) %>%
  inner_join(clean_vcf[c(1:3)], by = c("marker_name" = "ID"))
LG2drop <- chrom_map_df %>% group_by(LG_name) %>% summarise(num=n()) %>% filter(num<10)
drop_markers <- chrom_map_df %>% filter(LG_name %in% LG2drop$LG_name) %>% .[,"marker_name", drop=TRUE]
# mapRIL5_chrom_final <- 
mapRIL5_chrom_final <-  mstmap(jittermap(drop.markers(mapRIL5_chrom, drop_markers)), 
                               bychr = TRUE, trace = TRUE, 
                               dist.fun ="kosambi", p.value = 2)
# visualise
heatMap(mapRIL5_chrom_final, lmax = 30)
map_df <- imap_dfr(mapRIL5_chrom_final$geno, ~tibble(LG_name = .y, 
                                       marker_name = markernames(mapRIL5_chrom_final, .y))) %>%
  inner_join(clean_vcf[c(1:3)], by = c("marker_name" = "ID")) %>% 
  mutate(chrom_match = map2_lgl(.x=LG_name, .y=`#CHROM`,.f=grepl)) # check if LG matches Chrom
# summarise rates of matching LG:Chrom at each chromosome
map_match_rates <- map_df %>% group_by(LG_name) %>% summarise(match_rate=sum(chrom_match)/n())
plot.map(mapRIL5_chrom_final)

inds <- row.names(mapRIL5_chrom_final$geno[[1]]$data)
markernames(mapRIL5_chrom_final)
#### Save files for QTL analysis ####
# Export map
gmap_map <- imap_dfr(mapRIL5_chrom_final$geno, ~tibble(marker = markernames(mapRIL5_chrom_final, .y),
                                                       chr = .y, pos=.x$map)) %>%
  write_csv(sprintf("./data/qtl2_files/Lentil_GBS_%s_gmap.csv", stacks_name))
# Save physical map for R/qtl2
# Load marker map (extracted from gstacks.fa)
# stacks_dir <- "../Analysis/ref_stacks/stacks2_population_04_12_2017"
# # Load marker map (extracted from gstacks.fa)
# marker_map <- read_tsv(file.path(stacks_dir,"tag_chrom.map")) %>% 
#   dplyr::rename(CHROM_POS=POS) 
# # left_join(marker_map, c("#CHROM" = "LOC")) %>% 
# # mutate(ID=paste0("tag",`#CHROM`, POS, sep="_"), `#CHROM` = CHROM, 
# #        POS=ifelse(STRAND=="+", CHROM_POS+POS, CHROM_POS-POS)) %>% 
# # dplyr::select(-one_of(colnames(marker_map)))
# 
# marker_map <- read_tsv(file.path(stacks_dir, "tag_chrom.map")) %>%
#   dplyr::rename(CHROM_POS=POS)
# phys_map <- geno_map %>% 
#   separate(marker, into = c("tag", "LOC", "Phys_POS"), remove=FALSE, convert=TRUE) %>% 
#   inner_join(marker_map, .) %>% mutate(pos=Phys_POS) %>%
#   dplyr::select(one_of(colnames(geno_map))) %>% arrange(chr, pos)
# write_csv(phys_map, "data/qtl2_files/Lentil_GBS_pmap.csv")
# Write genotype table
genos <- as.data.frame(t(vcfR::extract.gt(vcfR::read.vcfR(clean_vcf_file)))) %>% 
  rownames_to_column %>% filter(rowname %in% inds) %>% .[c("rowname", markernames(mapRIL5_chrom_final))]
genos %>% filter(!grepl("RF$", rowname)) %>%
  write_csv(sprintf("./data/qtl2_files/Lentil_GBS_%s_genos.csv", stacks_name))
# Write founder table
genos %>% filter(grepl("RF$", rowname)) %>%
  write_csv(sprintf("./data/qtl2_files/Lentil_GBS_%s_founder_genos.csv", stacks_name))

# modify YAML files (need to have the templates in place)
yamls <- list.files("./data/qtl2_files/", "Lentil_GBS_LG.+\\d+dpi.yaml", full.names = TRUE)
for (y in yamls){
  
  system2("sed", args=c(glue::glue("-r 's/Lentil_GBS([A-z_]+.csv)/Lentil_GBS_{stacks_name}\\1/g'"),y),
          stdout = sub("LG", sprintf("%s", stacks_name), y))
}
  