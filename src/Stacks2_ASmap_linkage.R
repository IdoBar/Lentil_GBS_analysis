devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# Install and load required packages from GitHub
# .fixdevtools()
# remove.packages("onemap")
pacman::p_load_gh(c("augusto-garcia/onemap")) # "thierrygosselin/radiator", "augusto-garcia/onemap", "IdoBar/onemap@patch-3"
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
vcf_file <- "../Stack2_vcf_files/stacks_M1m4n2_populations.snps.vcf"  #
group_markers <- "CHR"
vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))
  
asma

#### ASMap ####
# cleaned file to read
clean_vcf_file <- "./data/intermediate_files/ref_stacks2_populations.snps.miss10.clean.vcf"
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

range_length <- 10
# p_value_range <- 5e-10*(5^(1:range_length))
p_value_range <- 1e-10*(2^(1:range_length))

# Create cluster with desired number of cores
registerDoFuture()
plan(multiprocess, workers=4)

asmap_thres_df <<- foreach(p = p_value_range, .combine = bind_rows) %dopar% {
    # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
    # l=3, rf=0.2
    # dummy <- length(onemap_file)
    testd <- R.utils::captureOutput(mstmap(raw_geno, pop.type = "RIL5", miss.thresh=0.1, p.value = p))
    # p.value = 5e-2
    
    data.frame(p.value = p, linkage_groups=as.numeric(stringr::str_extract(testd[1], "\\d+$")))
    
  }

asmap_thres_list <<- foreach(p = p_value_range, 
          .final = function(x) setNames(x, paste0("asmap_p", seq_along(p_value_range)))) %dopar% {
  # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
  # l=3, rf=0.2
  # dummy <- length(onemap_file)
  mstmap(raw_geno, pop.type = "RIL5", miss.thresh=0.1, p.value = p, bychr = FALSE) 
  # p.value = 5e-2 , objective.fun = "ML"
  
}
# Check length of groups
sapply(asmap_thres_list, function(l) qtl::chrlen(l))
# Check number of markers per groups
sapply(asmap_thres_list, function(l) qtl::nmar(l)) 

# selected p7 (from range)
selected_p <- 7
asmap <- asmap_thres_list[[selected_p]]
LogMsg(glue::glue("Selected p.value is: {p_value_range[selected_p]}"))
LogMsg(glue::glue("Number of markers in each linkage group: {qtl::nmar(asmap)}"))
LogMsg(glue::glue("Length of each linkage group: {qtl::chrlen(asmap)}"))




gc <- genClones(asmap, tol = 0.95)
if ("cgd" %in% names(gc)) {
  mapBC2 <- fixClones(asmap, gc$cgd, consensus = FALSE)
} else {mapBC2 <- asmap}
profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Remove highly distorted markers
mm <- statMark(mapBC2, stat.type = "marker")$marker$AA
LogMsg(glue::glue("Number of distorted markers:: {sum(mm > 0.7 | mm < 0.3)}"))
dm <- markernames(mapBC2)[(mm > 0.7) | (mm < 0.3)]
mapBC3 <- drop.markers(mapBC2, dm)
mapBC4 <- jittermap(mapBC3)
profileMark(mapBC4, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Pull markers that are co-located
mapDHs <- pullCross(mapBC4, type = "co.located")
# Pull markers that have significant segregation distortion with a p-value less than 0.02 and have a missing
# value proportion greater than 0.03.
# mapDHs <- pullCross(mapDHs, type = "seg.distortion", pars =
#                       list(seg.thresh = 0.02))
# mapDHs <- pullCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.08))
profileMark(mapDHs, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
nmar(mapDHs)
# Visualise data
heatMap(mapDHs, lmax = 20)
# Profile the statistics for the genotypes across
profileGen(mapDHs, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id =
             "Genotype", xo.lambda = 50, layout = c(1, 3), lty = 2)


mapBC5 <- mstmap(mapBC4, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = p_value_range[7])
profileMark(mapBC5, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")

qtl::nmar(mapBC5)
qtl::nmar(mapBC4)

