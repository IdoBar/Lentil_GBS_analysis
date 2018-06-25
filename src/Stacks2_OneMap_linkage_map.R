#!/bin/env Rscript
# devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# Install and load required packages from GitHub
# .fixdevtools()
# remove.packages("onemap")
options(unzip = "/usr/bin/unzip")
Sys.setenv(TAR = "/bin/tar")
# pacman::p_load_gh(c("IdoBar/onemap@patch-3")) # "thierrygosselin/radiator", "augusto-garcia/onemap"
# vignette("Inbred_Based_Populations") # OneMap tutorial
pacman::p_load_gh("augusto-garcia/onemap") # 
# devtools::install_github('hadley/ggplot2')
# pacman::p_load(Rhtslib)
# Install and load needed packages
package_list <- c("tidyverse", "qtl", "ASMap", "RColorBrewer", "Rsamtools", "doFuture", 
                  "foreach") # "radiator", "plotly"
pacman::p_load(char=package_list)

# Onemap requires R>=3.4 and Rhtslib and zlibbioc - install with conda
# git_packs <- c("augusto-garcia/onemap")
# install.deps(git_packs, repo = "git")
# devtools::install_github(git_packs)

options(stringsAsFactors=FALSE)
#### parse command line arguments ####
usage <- 'usage: Lentil_GBS_Stacks2_linkage_map_unix.R [options] <vcf_file>

options:
 -l, --loci <loci_rate_list>  String of rates of missing variant calls (in decreasing order) used for gradual filtration of the loci [default: 0.2,0.15,0.1].
 -g, --genotype_miss <geno_rate>   Rate of valid variant calls at each loci used to filter genotypes (samples) [default: 0.8].
 -n, --ncpus <ncores>              Number of cores to use for parallel processing [default: 8].
 -t, --tasks <tasks>               Up to which task to run (1. vcf_filtration, 2. find_LOD_rf, 3. group_LGs, 4. order_LG_markers) [default: 4].
 --lod <lod>                       LOD score to use in the final map construction [default: 6].
 --max_rf <rf>                     Recombinant Factor to use in the final map construction [default: 0.3].
 --map <map-type>                  Indicates the function that should be used, which can be "kosambi" or "haldane" [default: kosambi].
 --group <group-by>                Whether to group markers by Chromosome (chr) or by linkage group (lg) [default: chr]. 
 -h, --help                        Print this help message.'
#  -i, --input <vcf_file>            Stacks2-derived VCF file to process.
opts <- docopt::docopt(usage)
vcf_file <- opts$vcf_file #  args[1]
loci_miss_rates <- as.numeric(unlist(stringr::str_split(opts$loci, pattern = "[ ,_]+")))
geno_rate <- as.numeric(opts$genotype_miss)
ncpus <- as.numeric(opts$ncpus) # args[2]
LOD <- as.numeric(opts$lod)
RF <- as.numeric(opts$max_rf)
map_type <- opts$map
group_markers <- toupper(opts$group) # group_markers <- "LG"
# print command arguments
command_args <- unlist(opts[!grepl("^[-<]", names(opts))]) # %>% 
LogMsg(sprintf("Starting processing vcf file (%s).\nScript initiated with the following options:\n%s", vcf_file, 
               paste(paste(names(command_args), command_args, sep = ": "), collapse = ", ")))
# cat(paste(paste(names(command_args), command_args, sep = ": "), collapse = ", "))
# ncpus <- 8
tasks_to_run <- as.numeric(opts$tasks) # length(tasks)
# args <- commandArgs(trailingOnly = TRUE)
completed_tasks <- 0
# Load vcf file
# stacks_dir <- ""

# vcf_file <- "../Stack2_vcf_files/stacks_M1m4n2_populations.snps.vcf"  #
vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))

#### Load the most recent Rdata file ####

rdata_file_details <- file.info(list.files("./data/intermediate_files", full.names = TRUE,
                                           sprintf('%s.*[%s].*\\.RData', vcf_basename, group_markers))) %>% 
  rownames_to_column("filename") %>%
  mutate(mtime=as.POSIXct(mtime)) %>% arrange(desc(mtime)) 

if (nrow(rdata_file_details)>0) {
  LogMsg(sprintf("Data from a previous run was identified and loaded ('%s').", rdata_file_details$filename[1]))
  load(rdata_file_details$filename[1])
  last_completed_task <- stringr::str_extract(tasks[completed_tasks], "[^\\(]+")
  LogMsg(sprintf("Previous run terminated after task %d (%s).", completed_tasks, last_completed_task))
} 

# Processing functions
vcf_filtration <- function(vcf_file, miss_rates = seq(0.2, 0.1, -0.05), geno_rate=0.8) {
  LogMsg(sprintf("Filtering vcf file (%s)", vcf_file))
  vcf_data <- read_tsv(vcf_file, comment = "##") #%>% 
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
  # Write back vcf
  vcf_basename <<- tools::file_path_sans_ext(basename(vcf_file))
  clean_vcf_file <<- filedate(sprintf("%s.miss%d",vcf_basename, m*100), 
                             ".clean.vcf","data/intermediate_files", dateformat = FALSE)
                              # sub(".clean.sorted.vcf", "_miss40.clean.vcf",vcf_file, fixed = TRUE))
  
  # unlink(clean_vcf_file)
  # Copy the header into a new file
  system2("grep", args=c("'^##'", vcf_file), stdout = clean_vcf_file)
  if (file.exists(clean_vcf_file)) write_tsv(clean_vcf, clean_vcf_file, append = TRUE, col_names = TRUE)
  
  completed_tasks <<- 1
}
# Import to OneMap 
# Try to read straight from vcf
# clean_vcf_file <- "data/intermediate_files/ref_stacks2_populations.snps.miss10.clean.vcf"
find_LOD_rf <- function(clean_vcf_file, LOD_range = 3, rf_range = seq(0.1,0.5, 0.05), nworkers=8){
  clean_vcf_basename <<- sub(".clean.vcf", "", basename(clean_vcf_file))
  onemap_file <<- onemap_read_vcfR(vcfR::read.vcfR(clean_vcf_file), cross = "ri self",
                                  parent1 = "ILL6002_RF", 
                                  parent2 = "ILL7537_RF")
  
  # Check for marker segregation
  ril_test <- test_segregation(onemap_file) 
  Bonferroni_alpha(ril_test)
  pdf(filedate(sprintf("%s_seg_test", clean_vcf_basename), ".pdf", "plots/param_estimation"), 
      width = 8, height = 6)
  plot(ril_test)
  dev.off()
  good_markers <- select_segreg(ril_test, distorted = FALSE, numbers = TRUE)
  
  LOD_sug <- round(suggest_lod(onemap_file)) # estimate LOD threshold
  twopts_ril <<- rf_2pts(onemap_file, LOD=LOD_sug, max.rf = 0.35)
  # Assigning markers to linkage groups
  mark_all_ril <<- make_seq(twopts_ril, good_markers)
  
  if (length(LOD_range)==1) LOD_range <- (LOD_sug-LOD_range):(LOD_sug+LOD_range) # 4:10
  # Create cluster with desired number of cores
  registerDoFuture()
  plan(multiprocess, workers=nworkers)
  LogMsg(sprintf("Starting testing LOD and max.rf parameters.", paste(LOD_range, collapse = ","),
                 paste(rf_range, collapse = ",")))
  LGs_df <<- foreach(l = LOD_range, .combine = bind_rows) %:% 
    foreach(rf = rf_range, .combine = bind_rows) %dopar% {
      # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
      # l=3, rf=0.2
      dummy <- length(onemap_file)
      LG <- onemap::group(mark_all_ril, max.rf = rf, LOD = l)
      data.frame(LOD=l, max.rf=rf, as.data.frame(table(LG$groups)))
      
    }
  LogMsg(sprintf("Finished testing LOD (%s) and max.rf (%s) parameters.", 
                 paste(LOD_range, collapse = ","), paste(rf_range, collapse = ",")))
  # LGs_df <- bind_rows(LGs_df, df)
  # stopCluster(cl)
  
  
  colnames(LGs_df) <- c("LOD", "max.rf", "LG", "Freq")
  write_tsv(LGs_df, filedate(sprintf("%s_LOD_rf_params",clean_vcf_basename), ".tsv", "./data/intermediate_files", dateformat=FALSE))
  # Manualy check results:
  # LOD_range = 4:10
  # LGs_df <- read_tsv("data/intermediate_files/stacks_M1m4n2_populations.snps.miss10_LOD_rf_params_18_06_2018.tsv")
  LGs_sum <- LGs_df  %>% group_by(LOD, max.rf) %>% summarise(n=n()) %>% 
    mutate(LOD_cat=factor(LOD, levels=as.character(LOD_range)))
  y_max <- ceiling(max(LGs_sum$n)/5)*5  # find nearest x5 multiplier above the highest value
  # LGs_df %>% filter(LOD>7)
  p <- ggplot(LGs_sum , aes(x=max.rf, y=n, color=LOD_cat)) + geom_line(lwd=1) + 
    scale_color_brewer(palette = "Paired") + scale_y_continuous(name="Number of linkage groups", limits=c(0,y_max))
  p
  ggsave(filename = filedate(sprintf("%s_LOD_rf_est", clean_vcf_basename), ".pdf", "./plots/param_estimation", dateformat=FALSE), 
         width = 8, height = 6)
  completed_tasks <<- 2
  # save the environment
  save.image(filedate(sprintf("%s.%s",  clean_vcf_basename, group_markers), ".RData", 
                      "./data/intermediate_files", dateformat = FALSE))
}

group_LG <- function(nworkers=8, maptype = "kosambi", l=5, rf=0.35){
  set_map_fun(type = maptype) # haldane
  LG_markers <<- group(mark_all_ril, max.rf = rf, LOD = l)
  # Check markers distribution in Linkage groups
  table(LG$groups)
  completed_tasks <<- 3
  # save the environment
  save.image(filedate(sprintf("%s.%s",  clean_vcf_basename, group_markers), ".RData", 
                      "./data/intermediate_files", dateformat = FALSE))
  
}

 # 6
  #0.25
order_LG_markers <- function(LG, nworkers=8){
  # Create cluster with desired number of cores
  registerDoFuture()
  plan(multiprocess, workers=nworkers)
  
  # Order the markers in each of the linkage groups
  LogMsg(sprintf("Determining order of all (%d) linkage groups, please wait...", LG$n.groups))
  LG_ord <<- foreach(g = 1:LG$n.groups, .final = function(x) setNames(x, paste0("LG", 1:LG$n.groups))) %dopar% {
    # Order the markers in the current LG
    dummy <- length(onemap_file) + length(twopts_ril)
    # g=3
    LGg_f2 <- onemap::make_seq(LG, g)
    # Construct the linkage map, by automatic usage of the try algorithm:
    # LogMsg(sprintf("Determining order of LG %d, please wait", g))
    LGg_f2_ord <- onemap::order_seq(input.seq = LGg_f2, n.init = 5,
                                    subset.search = "twopt",
                                    twopt.alg = "rcd", THRES = 3,
                                    touchdown = TRUE)
    # LogMsg(sprintf("Finished ordering LG %d!", g))
    # Get the order with all markers:
    LGg_ordered <- onemap::make_seq(LGg_f2_ord, "force")
    pdf(filedate(sprintf("%s_%s_%s_rf_graph", clean_vcf_basename, group_markers, paste0("LG", g)), 
                 ".pdf", "./plots/temp_map_figures", dateformat=FALSE), 
        width = 20, height = 16)
    onemap::rf_graph_table(LGg_ordered, inter=FALSE)
    dev.off()
    return(LGg_ordered)
  }
  LogMsg(sprintf("Finished orderering all (%d) linkage groups.", LG$n.groups))
  # names(LG_list) <- paste0("LG", 1:LG$n.groups)
  completed_tasks <<- 4
  # save the environment
  save.image(filedate(sprintf("%s.%s",  clean_vcf_basename, group_markers), ".RData", 
                      "./data/intermediate_files", dateformat = FALSE))
}


group_CHR <- function(nworkers=8, maptype = "kosambi", l=5, rf=0.35){
  Chroms <- unique(gsub("_scaffold.+$", "", twopts_ril$CHROM))
  # Create cluster with desired number of cores
  registerDoFuture()
  plan(multiprocess, workers=nworkers)
  
  # Order the markers in each of the linkage groups
  LogMsg(sprintf("Grouping markers in each chromosome (%s-%s), please wait...", Chroms[1], Chroms[length(Chroms)]))
  Chr_list <<- foreach(c = Chroms, .final = function(x) setNames(x, Chroms)) %dopar% {
    # Make sure that the objects are available to the cluster workers
    dummy <- length(twopts_ril) + length(onemap_file)
    # Order the markers in the current Chromosome
    # c=LcChr3
    Chr_seq <- onemap::make_seq(twopts_ril, c)
  }
  CHR_markers <<- onemap::group_seq(twopts_ril, seqs=Chr_list, unlink.mks = mark_all_ril, rm.repeated = TRUE)
  LogMsg(sprintf("Finished grouping markers in all (%d) chromosomes.", length(Chroms)))
  # names(Chr_list) <- Chroms
  completed_tasks <<- 3
  # save the environment
  save.image(filedate(sprintf("%s.%s",  clean_vcf_basename, group_markers), ".RData", 
                      "./data/intermediate_files", dateformat = FALSE))
}

order_CHR_markers <- function(CHR_mks, nworkers=8){
  # CHR_mks <- CHR_markers
  Chroms <- names(CHR_mks$sequences)
  LogMsg(sprintf("Ordering markers in each chromosome (%s-%s), please wait...", Chroms[1], Chroms[length(Chroms)]))
  # Create cluster with desired number of cores
  registerDoFuture()
  plan(multiprocess, workers=nworkers)
  CHR_ord <<- foreach(c = Chroms, 
                      .final = function(x) setNames(x, Chroms)) %dopar% {
      dummy <- length(onemap_file) + length(twopts_ril) # + 
      CHR1_ord <- onemap::order_seq(CHR_mks$sequences[[c]])
      CHR1_frame <- onemap::make_seq(CHR1_ord, "force")
      pdf(filedate(sprintf("%s_%s_%s_rf_graph", clean_vcf_basename, group_markers, c), 
                   ".pdf", "./plots/temp_map_figures", dateformat=FALSE), 
          width = 20, height = 16)
      onemap::rf_graph_table(CHR1_frame, inter=FALSE)
      dev.off()
      return(CHR1_frame)
  }
  
  LogMsg(sprintf("Finished ordering markers in all (%d) chromosomes.", length(Chroms)))
  # names(Chr_list) <- Chroms
  completed_tasks <<- 4
  # save the environment
  save.image(filedate(sprintf("%s.%s",  clean_vcf_basename, group_markers), ".RData", 
                      "./data/intermediate_files", dateformat = FALSE))
}


# identify recent modified file and load it
#### Run tasks ####
# set tasks
tasks <- c("vcf_filtration(vcf_file, miss_rates = loci_miss_rates, geno_rate = geno_rate)", 
           "find_LOD_rf(clean_vcf_file, nworkers = ncpus)",
	sprintf("group_%s(nworkers = ncpus, maptype = map_type, l = LOD, rf = RF)", group_markers), 
	sprintf("order_%s_markers(%s_markers, nworkers = ncpus)", group_markers,group_markers))
# Chr_tasks <- c("vcf_filtration(vcf_file, miss_rates = loci_miss_rates, geno_rate = geno_rate)", "find_LOD_rf(clean_vcf_file, nworkers = ncpus)","group_LGs(mark_all_ril, maptype = map_type, l = LOD, rf = RF)", "order_CHR_markers(CHR_markers, nworkers = ncpus)") 

# run task 
while ((completed_tasks+1)<=tasks_to_run) {
  eval(parse(text=tasks[completed_tasks+1]))
}
  # see if get(text[completed_tasks+1]) works similarly

ordered_markers_list <- get(sprintf("%s_ord", group_markers))
# View linkage map
pdf(filedate(sprintf("%s_%s_linkage_map_draft", clean_vcf_basename, group_markers), ".pdf", "./plots"), 
    width = 10, height = 8)
draw_map(ordered_markers_list, names = TRUE, grid = FALSE, cex.mrk = 0.7)
dev.off()
# Number of markers in each Chromosome:
LogMsg(sprintf("Number of markers in each linkage group: %s", 
       paste(sapply(ordered_markers_list, function(lg) length(lg$seq.num)), collapse = ", ")))

# match chromosome positions of markers:
markers <- colnames(onemap_file$geno)
temp_map_file <- tempfile()
onemap::write_map(ordered_markers_list,temp_map_file)
temp_map <- read_delim(temp_map_file, delim = " ",
           col_names = c("LG", "marker", "pos")) %>%
  mutate(Phys_Chrom=onemap_file$CHROM[match(marker, markers)], 
         Phys_Pos=onemap_file$POS[match(marker, markers)], marker_num=match(marker, markers)) %>% 
  write_csv(filedate(sprintf("%s_%s_linkage_map_draft", clean_vcf_basename, group_markers), 
                     ".csv", "./data/intermediate_files"))


# Manually adjust markers at edges of LGs:
# View the LOD scores and RF of LG4
# rf_graph_table(LG_list$LG4)
# LG_list$LG4 <- drop_marker(LG_list$LG4, c(960:967))
# LG4_ord <- order_seq(input.seq = LG_list$LG4, n.init = 5,
#                         subset.search = "twopt",
#                         twopt.alg = "rcd", THRES = 3,
#                         touchdown = TRUE)
# LogMsg(sprintf("Finished ordering LG %d!", g))
# # Get the order with all markers:
# LG4_seq <- make_seq(LG4_ord, "force")
# # Checking for better orders: (seems like the mapping order is pretty good)
# # ripple_seq(LG4_seq, ws = 5)
# # # try to add them again
# # LG4_add <- try_seq(LG4_seq, 967)
# # Insert back into LG_list
# LG_list$LG4 <- LG4_seq
# 
# # View the LOD scores and RF of LG2
# rf_graph_table(LG_list$LG2)
# LG_2_trim <- drop_marker(LG_list$LG2, c(438:439, 531:534))
# LG2_ord <- order_seq(input.seq = LG_2_trim, n.init = 5,
#                      subset.search = "twopt",
#                      twopt.alg = "rcd", THRES = 3,
#                      touchdown = TRUE)
# LogMsg(sprintf("Finished ordering LG %d!", g))
# # Get the order with all markers:
# LG2_seq <- make_seq(LG2_ord, "force")
# # Insert back into LG_list
# LG_list$LG2 <- LG2_seq
# # View linkage map
# draw_map(LG_list, names = FALSE, grid = FALSE, cex.mrk = 0.7)
# 
# #### Save files for QTL analysis ####
# # Export linkage map
# write_map(LG_list,filedate("Lentil_GBS_linkage_map_order", ".map", "data", dateformat = FALSE))
# 
