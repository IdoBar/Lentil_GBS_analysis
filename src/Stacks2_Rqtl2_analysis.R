devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# Install R/qtl2
if (!require("qtl2")) install.packages("qtl2", repos="https://rqtl.org/qtl2cran")

# install.packages("mixtools")
# devtools::install_github("KonradZych/phenotypes2genotypes")

# Install and load needed packages
package_list <- c("tidyverse", "qtl", "RColorBrewer", "qtl2") # "radiator", "qtl2"
pacman::p_load(char = package_list)



options(stringsAsFactors = FALSE)


#### Analysis in R/qtl2  ####
stacks_name <- "ref_stacks2"
yaml_files <- dir("data/qtl2_files", "Lentil_GBS_.+\\.yaml", full.names = TRUE)
Chr_search="LG"
yaml_files <- yaml_files[grepl(stacks_name, yaml_files)]
time_points <- as.character((1:4)*7)
thresh_df <- NULL
LOD_scores <- vector(mode = "list", length = length(time_points))
peaks_df <- NULL
qtl_list <- vector(mode = "list", length = length(time_points))
for (dpi in time_points){
  # dpi=time_points[1]
  yml <-yaml_files[grepl(dpi, yaml_files)]
  LogMsg(sprintf("Processing YAML file: %s", yml))
  gbs_data <- read_cross2(yml)
  # insert pseudomarkers into the genetic map
  
  if (!exists("gen_map")) {
    # gen_map=gbs_data$gmap
    gen_map <- insert_pseudomarkers(gbs_data$gmap, step=1)
    gen_map <- reduce_map_gaps(gen_map)
    gen_map_df <- lapply(names(gen_map), function(chr) data.frame(Locus=names(gen_map[[chr]]), 
                                                              Position=gen_map[[chr]]) %>%  
                       mutate(Group=chr)) %>% map_df(bind_rows)
  }
  
  # calculate the QTL genotype probabilities
  # gen_map=gbs_data$gmap
  pr <- calc_genoprob(gbs_data, gen_map, error_prob=0.002, cores=2)
  # Calculate error LOD 
  # calc_errorlod()
  # Calculate kinship matrix (needed for LMM)
  kinship_loco <- calc_kinship(pr, "loco", cores=2)
  # perform a genome scan by LMM analysis  
  qtl_list[[dpi]] <- scan1(pr, gbs_data$pheno, kinship_loco, cores=2)
  LOD_scores[[dpi]] <- as.data.frame(qtl_list[[dpi]]) %>% rownames_to_column("Locus") %>% # [[dpi]]
    gather(key="Trait", value="LOD", Stem_breakage_percent:Leaf_lesion_percent) %>% 
    left_join(., gen_map_df) %>% dplyr::select(Position, Group, LOD, Trait) 
  # calculate LOD threshold based on permutation step
  # perm_lod <- scan1perm(pr, gbs_data$pheno, kinship_loco, n_perm = 1000,cores=2)
  # thres_sum <- data.frame(summary(perm_lod)) %>% rownames_to_column("alpha") %>% mutate(dpi=dpi)
  # thresh_df <- bind_rows(thresh_df, thres_sum)
  
  # Find peaks
  peaks <- find_peaks(qtl_list[[dpi]], gen_map, threshold=3, drop=1, peakdrop=2) %>% mutate(dpi=dpi)
  peaks_df <- bind_rows(peaks_df, peaks) # thres_sum$Leaf_lesion_percent
  # get physical position at location - need to add physical map (pmap) to the cross - check R/qtl
  # peak_Mbp <- max(qtl_list[[dpi]], gbs_data$pmap)$pos
  
}


##### export to MapChart ####
# Save all information to an mct file (for MapChart)
mc_file <- filedate("Lentil_ascochyta_QTL_map", ".mct", "QTL_results")
unlink(mc_file)
raw_map_df <- names(gbs_data$gmap) %>%
  map(function(chr) data.frame(Locus=names(gbs_data$gmap[[chr]]),
                               Position=gbs_data$gmap[[chr]]) %>%
        mutate(Group=paste0("Chr", chr))) %>% map_df(bind_rows)
peaks_df <- unfactor(peaks_df)
LGs <- unique(gen_map_df$Group)

traits <- unique(peaks_df$lodcolumn)
# names(time_points) <- c(0:2,7)
time_fill <- setNames(c(0:2,7), unique(peaks_df$dpi))
for (i in seq_along(LGs)){
  # i=1
  ch <- LGs[i]
  # Write group name to file
  write(paste("chrom", ch), mc_file, append = file.exists(mc_file))
  write_tsv(gen_map_df %>% filter(Group==ch) %>% select(-Group), mc_file, append = TRUE, col_names = FALSE)
  if (peaks_df %>% filter(chr==i) %>% nrow(.) > 0)  write(sprintf("\nqtls"),  mc_file, append =TRUE)
  for (t in seq_along(traits)){
    for (d in names(time_fill)){
    
      # d=1
      # t=1
      trait <- traits[t]
      # save LOD scores to files
      lod_intervals <- peaks_df %>% filter(chr==i, lodcolumn==trait, dpi==d) %>% select(ci_lo, ci_hi) #%>%
        # as.character()
      if (nrow(lod_intervals)>0) {
        interval_string <- paste(rep(round(lod_intervals, digits=3), each=2), collapse = " ")
        # Use filenames in MCT file
        write(sprintf("%s_%s %s C%d F%d", trait, d, interval_string, t, time_fill[d]),  mc_file, append =TRUE) # C4 F5
        # write(sprintf("%s_%s QTL_results/%s_%sdpi_%s.lod C1 L1 S0\n", t, dpi, t, dpi, ch),  
      }
      
      #       mc_file, append =TRUE) 
    }
    #   
  }
  write(sprintf("\n"),  mc_file, append =TRUE)
}

xlsx::write.xlsx(as.data.frame(peaks_df), 
                 filedate(glue::glue("LOD_peaks_{stacks_name}_mapping"), 
                          ".xlsx", "QTL_results"), row.names = FALSE)
# xlsx::write.xlsx(as.data.frame(peaks_df), filedate("LOD_peaks_Chr_mapping", ".xlsx", "data"), row.names = FALSE)
# thresh_df
peaks_df %>% group_by(chr) %>% summarise(max_lod=max(lod))

# Compare the LOD scores for each trait with each method
# color <- c("slateblue", t_col("violetred", percent = 0))
colors <- adjustcolor(RColorBrewer::brewer.pal(4, "Set1"), alpha.f = 0.65)

#par(mar=c(4.1, 4.1, 1.6, 1.1))
par(mfrow=c(2,2))
# lod_vars <- paste0("out_pg_", time_points)
ymx <- max(peaks_df$lod)

pdf(filedate(glue::glue("Lentil_GBS_LOD_{stacks_name}"), ".pdf", "plots"), width = 8, height = 6)

for(i in 1:4) {
  # i=4
  plot(qtl_list$`7`, gen_map, lodcolumn=i, col=colors[1], main=colnames(gbs_data$pheno)[i], 
       xlab="Linkage Group",
       ylim=c(0, ymx*1.02))
  abline(h=thresh_df[1,i+1], col=colors[1], lty="dashed", lwd=0.8)
  for (o in 2:length(time_points)){
    # o=3
    
    plot(qtl_list[[time_points[o]]], gen_map, lodcolumn=i, col=colors[o],
                            add = TRUE)
    abline(h=thresh_df[o,i+1], col=colors[o], lty="dashed", lwd=0.8)
  }
  
  legend("topright", lwd=2, col=colors, paste(time_points, "dpi"), bg="gray90") # "H-K", , lty=c(1)
}
dev.off()


# Plot each trait and time separately
for(i in 1:4) {
  for (o in 1:length(time_points)){
  # i=4
    plot(qtl_list[[time_points[o]]], gen_map, lodcolumn=i, col=colors[o], main=colnames(gbs_data$pheno)[i],
                            xlab='Linkage Group',ylim=c(0, ymx*1.02))
    abline(h=thresh_df[o,i+1], col=colors[o], lty="dashed", lwd=0.8)
   legend("topright", lwd=2, col=colors, paste(time_points, "dpi"), bg="gray90")
    # o=3
    
  }
  
   # "H-K", , lty=c(1)
}
# Focus on one chromosome
chrom <- "LG3"
pdf(filedate(sprintf("Lentil_GBS_LOD_%s", chrom), ".pdf", "plots"), width = 8, height = 6)
for(i in 1:4) {
  plot(qtl_list[[1]], map, lodcolumn=i, col=colors[1], chr=chrom, main=colnames(gbs_data$pheno)[i], xlab=sprintf("Linkage Group %s", chrom),
       ylim=c(0, ymx*1.02))
  abline(h=thresh_df[1,i+1], col=colors[1], lty="dashed", lwd=0.8)
  for (o in 2:length(lod_vars)){
    plot(qtl_list[[time_points[o]]], map, lodcolumn=i, col=colors[o],chr=chrom,add = TRUE)
    abline(h=thresh_df[o,i+1], col=colors[o], lty="dashed", lwd=0.8)
  }
  
  legend("topleft", lwd=2, col=colors, paste(time_points, "dpi"), bg="gray90") # "H-K", , lty=c(1)
}
dev.off()

#### Estimating QTL effect ####

# Read the data for 14 dpi
gbs_data <- read_cross2(yaml_files[grepl("14", yaml_files)])
pr <- calc_genoprob(gbs_data, gen_map, error_prob=0.002, cores=2)
# Calculate kinship matrix (needed for LMM)
kinship_loco <- calc_kinship(pr, "loco", cores=2)
# get the estimated effects on chromosome 2 for the liver phenotype, we'd do the following:
LG3eff <- scan1coef(pr[,"3"], gbs_data$pheno[,"Leaf_lesion_percent"], kinship = kinship_loco[["3"]])
# Plot the SNP effect across the linkage group
col <- brewer.pal(2, "Set1")#c("slateblue", "green3", "violetred")
plot(LG3eff, gen_map["3"], columns=1:2, col=col)
last_coef <- unclass(LG3eff)[nrow(LG3eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

# Find markers under the QTL region
region <- peaks_df %>% filter(lod==max(lod))
lg_map <- gen_map[[region$chr]]
sel_markers <- lg_map[lg_map>=region$ci_lo & lg_map<=region$ci_hi]
orig_markers <- sel_markers[!grepl(".loc\\d+", names(sel_markers))]
# get inferred genotypes at the LG3 QTL
pdf(filedate("Genotype-phenotype_dist_LcChr2_QTL_horz", ".pdf", "plots"), width=10, height=10)
par(mfrow=c(3,2))
for (m in orig_markers){
  g <- maxmarg(pr, gen_map, chr=3, pos=m, return_char=TRUE)
  # plot phenotrype against genotype
  plot_pxg(g, gbs_data$pheno[,"Leaf_lesion_percent"], ylab="Leaf Lesion Coverage", SEmult=2, swap_axes=TRUE)
  title(main=sprintf("Genotype distribution at QTL3 (%.2fcM)", m))
}
dev.off()

# Find other markers at that range (from the genome mapping)
# Translate linkage group to corresponding chromosome
lg2chr <- set_names(c("LcChr1", "LcChr3", "LcChr2", "LcChr6", "LcChr4", "LcChr3", "LcChr4", 
                      "LcChr7", "LcChr5", "LcChr7"), names(gen_map))
# extract the significant markers from the genomic map
# Load marker map (extracted from gstacks.fa)
marker_map <- read_tsv("data/tag_chrom.map") %>% 
  dplyr::rename(CHROM_POS=POS) 

file_prefix <- "Lentil_GBS"


#### Import to R/qtl  ####

# In order to import to R/qtl, the cross type in the raw file should be changed to 'bc'
dat1 <- read.cross("mm", dir="data", 
                   file="gstacks_no_RL_NS_GQ_DP_Samples_filtered.rqtl.pheno.21dpi.raw", 
                   mapfile = "Lentil_GBS_linkage_chromosome_map_22_01_2018.map") # Lentil_GBS_linkage_groups_map_22_01_2018.map
summary(dat1)
# Then convert it to RIL
ril_dat <- convert2riself(dat1)
summary(ril_dat)
plot(ril_dat)
#### Exploratory analysis of the map ####
plotMissing(ril_dat)
# par(mfrow=c(1,2), las=1)
# identify missing information
plot(ntyped(ril_dat), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(ril_dat, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
# Identify duplicated genomes
cg <- comparegeno(ril_dat)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
g <- pull.geno(ril_dat)
table(g[38,], g[39,])
# Identify and remove duplicated markers
print(dup <- findDupMarkers(ril_dat, exact.only=FALSE))
ril_dat <- drop.markers(ril_dat, unlist(dup))
# Look for distorted markers
gt <- geno.table(ril_dat)
gt[gt$P.value < 0.05/totmar(ril_dat),] # None
# Asses RF
mapthis <- est.rf(ril_dat)
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6) # 0.35, 6
table(lg[,2])
plotRF(lg, alternate.chrid=TRUE)
# Order markers
mapthis <- orderMarkers(mapthis, chr=4)
rip4 <- ripple(mapthis, chr=4, window=6)
# Check likelihood (look for higher LOD)
rip4lik <- ripple(mapthis, chr=4, window=4, method="likelihood",
                                     error.prob=0.005)
# Check error LOD scores
ril_dat <- calc.errorlod(ril_dat)
# Remove linkage groups 8, 10:14 (contain very few markers) - for LG analysis
ril_dat <- pheno2geno::removeChromosomes(ril_dat, as.character(c(8,10:14)), verbose=TRUE)
# Check map
pdf("plots/Lentil_GBS_linkage_chromosome_map_22_01_2018.pdf", width = 8, height = 6)
plotMap(pull.map(ril_dat)) # ,est.map(ril_dat))
dev.off()
# plotMap(est.map(ril_dat))
summary(pull.map(ril_dat))
plotRF(ril_dat)





