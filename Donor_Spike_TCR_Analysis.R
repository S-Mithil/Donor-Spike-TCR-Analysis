# =============================================================================
# Donor_Spike_TCR_Analysis.R | 
# =============================================================================
# Data: https://clients.adaptivebiotech.com/pub/soni-2026-s
# =============================================================================

required_packages <- c("readr","readxl","dplyr","tidyr","ggplot2",
                       "ggalluvial","RColorBrewer","openxlsx","scales")
missing_packages <- setdiff(required_packages,
  required_packages[vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)])
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages,collapse=", "))
suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr)
  library(ggplot2); library(ggalluvial); library(RColorBrewer)
  library(openxlsx); library(scales)
})
project_dir <- Sys.getenv("TCMR_PROJECT_DIR", unset=getwd())
raw_dir <- file.path(project_dir,"data","raw")
processed_dir <- file.path(project_dir,"data","processed")
results_dir <- file.path(project_dir,"results")
figures_dir <- file.path(results_dir,"figures")
tables_dir <- file.path(results_dir,"tables")
for (p in c(processed_dir,results_dir,figures_dir,tables_dir))
  dir.create(p,recursive=TRUE,showWarnings=FALSE)
repertoire_file <- file.path(raw_dir,"Pt1006_All_Rearrangements_Ambiguous_removed_plus_30M_54M.csv")
pre_cd4_file <- file.path(raw_dir,"CD4_Direct_PreDRTCCs.csv")
pre_cd8_file <- file.path(raw_dir,"CD8_Direct_PreDRTCCs.csv")
spike_file <- file.path(raw_dir,"Pt1006_Spike_Reactive_Clones_30M_33-7M_54M.xlsx")
for (p in c(repertoire_file,pre_cd4_file,pre_cd8_file,spike_file))
  if (!file.exists(p)) stop("Required input missing: ",p)

min_functional_frequency <- 1e-5
fold_enrichment_cutoff <- 2
comparator_fold_expansion <- 5
comparator_min_templates <- 8
comparator_min_frequency <- 1e-5
# Verify these sample-specific productive template totals against immunoSEQ.
cd8_30m_blood_total_templates <- 43614
bulk_fcb_total_templates <- 13765
alluvial_top_n_each_timepoint <- 20
alluvial_mean_percent_cutoff <- 0.75
representative_crtccs <- c("CASSLGRLRSPLHF","CASWSNRGQETQYF",
  "CASSEFAPRGSGDYEQYF","CASSQSTGDWDTQYF","CASSSGGLEGEQYV")
special_clone_colors <- c("CASSLGRLRSPLHF"="#E31A1C",
  "CASWSNRGQETQYF"="#1F78B4","CASSSGGLEGEQYV"="#6A3D9A")

# --- Helpers -----------------------------------------------------------------
clean_cdr3 <- function(x) { x <- trimws(as.character(x)); x[x==""] <- NA_character_; x }
clean_set <- function(x) unique(x[!is.na(x) & nzchar(x)])
assert_columns <- function(d,cols,label="data") {
  miss <- setdiff(cols,names(d))
  if (length(miss)) stop(label," missing columns: ",paste(miss,collapse=", "))
}
read_cdr3_csv <- function(path) {
  x <- readr::read_csv(path,show_col_types=FALSE)
  if (!"CDR3" %in% names(x)) names(x)[1] <- "CDR3"
  clean_set(clean_cdr3(x$CDR3))
}
truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x==1)
  toupper(trimws(as.character(x))) %in% c("TRUE","T","1","YES","Y")
}
make_post_drtcc <- function(d,stim,base) {
  assert_columns(d,c("CDR3",stim,base),"repertoire")
  clean_set(d$CDR3[d[[stim]] >= min_functional_frequency &
    (d[[base]]==0 | d[[stim]] > fold_enrichment_cutoff*d[[base]])])
}
write_table <- function(x,name) readr::write_csv(x,file.path(tables_dir,name),na="")
competition_rank <- function(x) {
  out <- rep(NA_integer_,length(x)); ok <- !is.na(x) & x>0
  out[ok] <- rank(-x[ok],ties.method="min"); out
}

# --- Section 1: cleaned repertoire matrix ------------------------------------
rep_raw <- readr::read_csv(repertoire_file,show_col_types=FALSE)
if (!"CDR3" %in% names(rep_raw)) names(rep_raw)[1] <- "CDR3"
rep_raw$CDR3 <- clean_cdr3(rep_raw$CDR3)
freq_cols <- grep("TCRB$",names(rep_raw),value=TRUE)
if (!length(freq_cols)) stop("No TCRB frequency columns")
rep_raw <- rep_raw |> filter(!is.na(CDR3),CDR3!="CDR3")
rep_raw[freq_cols] <- lapply(rep_raw[freq_cols],function(x) {
  y <- suppressWarnings(as.numeric(x)); y[is.na(y)] <- 0; y
})
repertoire <- rep_raw |> group_by(CDR3) |>
  summarise(across(all_of(freq_cols),~sum(.x,na.rm=TRUE)),.groups="drop")
readr::write_csv(repertoire,file.path(processed_dir,"repertoire_matrix_clean.csv"))

# --- Section 2: functional clone definitions ---------------------------------
pre_cd4 <- read_cdr3_csv(pre_cd4_file)
pre_cd8 <- read_cdr3_csv(pre_cd8_file)
missing_pre_cd4 <- setdiff(pre_cd4,repertoire$CDR3)
missing_pre_cd8 <- setdiff(pre_cd8,repertoire$CDR3)
write_table(data.frame(Lineage=c(rep("CD4",length(missing_pre_cd4)),
  rep("CD8",length(missing_pre_cd8))),CDR3=c(missing_pre_cd4,missing_pre_cd8)),
  "Audit_missing_pre_DRTCC.csv")
if (length(missing_pre_cd4)||length(missing_pre_cd8))
  stop("Pre-transplant DRTCC identifiers missing from cleaned matrix; inspect audit file")

drtcc <- list(
 CD4=list(PreTx=pre_cd4,
   `30M`=make_post_drtcc(repertoire,"30M_CFSE_Low_CD4_TCRB","30M_Unstim_CD4_TCRB"),
   `33.7M`=make_post_drtcc(repertoire,"33-7M_CFSE_low_CD4_TCRB","30M_Unstim_CD4_TCRB"),
   `54M`=make_post_drtcc(repertoire,"54M_CFSE_Low_CD4_TCRB","54M_Unstim_CD4_TCRB")),
 CD8=list(PreTx=pre_cd8,
   `30M`=make_post_drtcc(repertoire,"30M_CFSE_Low_CD8_TCRB","30M_Unstim_CD8_TCRB"),
   `33.7M`=make_post_drtcc(repertoire,"33-7M_CFSE_low_CD8_TCRB","30M_Unstim_CD8_TCRB"),
   `54M`=make_post_drtcc(repertoire,"54M_CFSE_Low_CD8_TCRB","54M_Unstim_CD8_TCRB")))
spike_cd4 <- readxl::read_excel(spike_file,sheet="CD4_spike_reactive")
spike_cd8 <- readxl::read_excel(spike_file,sheet="CD8_spike_reactive")
spike_cd4$CDR3 <- clean_cdr3(spike_cd4$CDR3)
spike_cd8$CDR3 <- clean_cdr3(spike_cd8$CDR3)
get_spike <- function(d,col) {
  assert_columns(d,c("CDR3",col),"curated spike sheet")
  clean_set(d$CDR3[truthy(d[[col]])])
}
srtcc <- list(
 CD4=list(`30M`=get_spike(spike_cd4,"CD4_hit_30M_spike"),
   `33.7M`=get_spike(spike_cd4,"CD4_hit_33_7M_spike"),
   `54M`=get_spike(spike_cd4,"CD4_hit_54M_spike")),
 CD8=list(`30M`=get_spike(spike_cd8,"CD8_hit_30M_spike"),
   `33.7M`=get_spike(spike_cd8,"CD8_hit_33_7M_spike"),
   `54M`=get_spike(spike_cd8,"CD8_hit_54M_spike")))
saveRDS(list(DRTCC=drtcc,SRTCC=srtcc),file.path(processed_dir,"functional_clone_sets.rds"))

# --- Section 3: Fig. 2B/D pairwise overlap and 2C/E frequencies --------------
mlr_cols <- list(PreTx=c(CD4="UPITT_1006_CFSElow_CD4_TCRB",CD8="UPITT_1006_CFSElow_CD8_TCRB"),
 `30M`=c(CD4="30M_CFSE_Low_CD4_TCRB",CD8="30M_CFSE_Low_CD8_TCRB"),
 `33.7M`=c(CD4="33-7M_CFSE_low_CD4_TCRB",CD8="33-7M_CFSE_low_CD8_TCRB"),
 `54M`=c(CD4="54M_CFSE_Low_CD4_TCRB",CD8="54M_CFSE_Low_CD8_TCRB"))
spike_cols <- list(`30M`=c(CD4="30M_SCoV2_S_CD4_TCRB",CD8="30M_SCoV2_S_CD8_TCRB"),
 `33.7M`=c(CD4="33-7M_SCoV2_S_CD4_TCRB",CD8="33-7M_SCoV2_S_CD8_TCRB"),
 `54M`=c(CD4="54M_SCoV2_S_CD4_TCRB",CD8="54M_SCoV2_S_CD8_TCRB"))
pair_summary <- list(); pair_clones <- list(); idx <- 0L
for (lin in c("CD4","CD8")) for (m in names(mlr_cols)) for (s in names(spike_cols)) {
  idx <- idx+1L
  ids <- intersect(drtcc[[lin]][[m]],srtcc[[lin]][[s]])
  mc <- unname(mlr_cols[[m]][lin]); sc <- unname(spike_cols[[s]][lin])
  sub <- repertoire |> filter(CDR3 %in% ids) |>
    transmute(CDR3,MLR_frequency_percent=100*.data[[mc]],
              Spike_frequency_percent=100*.data[[sc]])
  pair_clones[[idx]] <- data.frame(Lineage=lin,MLR=m,Spike=s,sub)
  pair_summary[[idx]] <- data.frame(Lineage=lin,MLR=m,Spike=s,
    MLR_DRTCC_count=length(drtcc[[lin]][[m]]),
    Spike_SRTCC_count=length(srtcc[[lin]][[s]]),Shared_count=length(ids),
    Cumulative_MLR_frequency_percent=sum(sub$MLR_frequency_percent),
    Cumulative_Spike_frequency_percent=sum(sub$Spike_frequency_percent))
}
pair_summary <- bind_rows(pair_summary); pair_clones <- bind_rows(pair_clones)
write_table(pair_summary,"Fig2BD_pairwise_overlap_summary.csv")
write_table(pair_clones,"Fig2BDE_pairwise_shared_clone_frequencies.csv")
write_table(bind_rows(
  pair_summary |> transmute(Lineage,MLR,Spike,Repertoire="MLR",
    Cumulative_productive_frequency_percent=Cumulative_MLR_frequency_percent),
  pair_summary |> transmute(Lineage,MLR,Spike,Repertoire="Spike",
    Cumulative_productive_frequency_percent=Cumulative_Spike_frequency_percent)),
  "Fig2CE_cumulative_shared_frequency.csv")
cd8_shared <- sort(unique(pair_clones$CDR3[pair_clones$Lineage=="CD8"]))

# --- Section 4: Fig. 2F alluvial (selection independent of reactivity) --------
blood_cols <- c(PreTx="UPITT_1006_UNSTIM_CD8_TCRB",`6M`="UPITT_1006_6mth_CD8_TCRB",
 `12M`="UPITT_1006_12mth_CD8_TCRB",`30M`="30M_Unstim_CD8_TCRB",
 `54M`="54M_Unstim_CD8_TCRB")
assert_columns(repertoire,unname(blood_cols))
top_union <- unique(unlist(lapply(blood_cols,function(col)
  repertoire |> arrange(desc(.data[[col]]),CDR3) |>
  slice_head(n=alluvial_top_n_each_timepoint) |> pull(CDR3))))
fig2f_wide <- repertoire |> filter(CDR3 %in% top_union) |>
  transmute(CDR3,PreTx=.data[[blood_cols["PreTx"]]],`6M`=.data[[blood_cols["6M"]]],
    `12M`=.data[[blood_cols["12M"]]],`30M`=.data[[blood_cols["30M"]]],
    `54M`=.data[[blood_cols["54M"]]]) |>
  mutate(Mean_percent=100*rowMeans(across(c(PreTx,`6M`,`12M`,`30M`,`54M`))),
    Max_frequency=pmax(PreTx,`6M`,`12M`,`30M`,`54M`),
    Shared_DRTCC_SRTCC=CDR3 %in% cd8_shared) |>
  filter(Mean_percent>=alluvial_mean_percent_cutoff) |>
  arrange(desc(Max_frequency),CDR3)
if(nrow(fig2f_wide)!=19) warning("Expected 19 alluvial clones; observed ",nrow(fig2f_wide))
clone_order <- fig2f_wide$CDR3
palette <- unique(c(RColorBrewer::brewer.pal(12,"Paired"),
  RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(8,"Set2"),
  RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(9,"Set1")))
palette <- palette[!toupper(palette) %in% toupper(unname(special_clone_colors))]
ordinary <- setdiff(clone_order,names(special_clone_colors))
clone_colors <- setNames(rep(NA_character_,length(clone_order)),clone_order)
clone_colors[ordinary] <- palette[seq_along(ordinary)]
for (id in intersect(clone_order,names(special_clone_colors))) clone_colors[id] <- special_clone_colors[id]
fig2f_long <- fig2f_wide |> select(CDR3,PreTx,`6M`,`12M`,`30M`,`54M`) |>
  pivot_longer(-CDR3,names_to="Timepoint",values_to="Frequency") |>
  mutate(Timepoint=factor(Timepoint,levels=names(blood_cols)),
    CDR3=factor(CDR3,levels=rev(clone_order)))
p_f <- ggplot(fig2f_long,aes(x=Timepoint,stratum=CDR3,alluvium=CDR3,
  y=Frequency,fill=CDR3)) +
  geom_alluvium(width=.25,knot.pos=.5,curve_type="xspline",
    colour="black",linewidth=.25) +
  geom_stratum(width=.25,colour="black",linewidth=.25,alpha=.35) +
  scale_fill_manual(values=clone_colors,guide="none") +
  scale_y_continuous(labels=scales::label_percent(accuracy=1),expand=c(0,0)) +
  labs(x=NULL,y="Cumulative productive frequency (%)") + theme_classic(base_size=11)
ggsave(file.path(figures_dir,"Fig2F_dominant_CD8_alluvial.pdf"),p_f,width=8.5,height=4.6)
ggsave(file.path(figures_dir,"Fig2F_dominant_CD8_alluvial.png"),p_f,width=8.5,height=4.6,dpi=600)
write_table(mutate(fig2f_wide,Plot_color=clone_colors[CDR3]),"Fig2F_alluvial_source_data.csv")

# --- Section 5: Fig. 2H validated mutually exclusive categories --------------
pre_drtcc <- drtcc$CD8$PreTx
post_drtcc_complete <- unique(unlist(drtcc$CD8[c("30M","33.7M","54M")],use.names=FALSE))
all_srtcc <- unique(unlist(srtcc$CD8,use.names=FALSE))
all_drtcc <- union(pre_drtcc,post_drtcc_complete)
pre_crtcc <- intersect(pre_drtcc,all_srtcc)
post_crtcc <- setdiff(intersect(post_drtcc_complete,all_srtcc),pre_crtcc)
all_crtcc <- union(pre_crtcc,post_crtcc)
# Pre-assignment takes precedence for donor-reactive-only clones as well.
pre_drtcc_only <- setdiff(pre_drtcc,all_srtcc)
post_drtcc_only <- setdiff(post_drtcc_complete,union(pre_drtcc,all_srtcc))
srtcc_only <- setdiff(all_srtcc,all_drtcc)
if (length(pre_crtcc)!=10 || length(post_crtcc)!=7 || length(all_crtcc)!=17)
  stop("CRTCC audit failed: expected 10 Pre / 7 Post / 17 unique; found ",
       length(pre_crtcc)," / ",length(post_crtcc)," / ",length(all_crtcc))
if (!setequal(all_crtcc,cd8_shared))
  stop("Fig. 2H CRTCC set differs from union of 12 pairwise CD8 comparisons")
# Comparator: 12M -> 30M after primary series, BEFORE booster.
# Uses 30M blood productive template total (not biopsy total).
rep_comparator <- repertoire |> mutate(
  Blood12=.data[["UPITT_1006_12mth_CD8_TCRB"]],
  Blood30=.data[["30M_Unstim_CD8_TCRB"]],
  Blood30_templates=round(Blood30*cd8_30m_blood_total_templates),
  Fold_12_to_30=case_when(Blood12>0~Blood30/Blood12,
    Blood12==0 & Blood30>0~Inf,TRUE~0))
comparator <- rep_comparator |> filter(Blood30>comparator_min_frequency,
  Blood30_templates>=comparator_min_templates,
  Fold_12_to_30>=comparator_fold_expansion,
  !CDR3 %in% union(all_drtcc,all_srtcc)) |> pull(CDR3)
fig2h_categories <- list(`Pre-DRTCC only`=pre_drtcc_only,
 `Post-DRTCC only`=post_drtcc_only,`SRTCC only`=srtcc_only,
 `Pre-CRTCC`=pre_crtcc,`Post-CRTCC`=post_crtcc,
 `nDRTCC-nSpike`=comparator)
cat_pairs <- combn(names(fig2h_categories),2,simplify=FALSE)
overlap_audit <- bind_rows(lapply(cat_pairs,function(p) data.frame(
  Category_1=p[1],Category_2=p[2],Overlap_count=length(intersect(
    fig2h_categories[[p[1]]],fig2h_categories[[p[2]]])))))
write_table(overlap_audit,"Audit_Fig2H_category_overlap.csv")
if (any(overlap_audit$Overlap_count!=0)) stop("Fig. 2H categories overlap; see audit CSV")
fcb_col <- "Upitt_1006_FcB9M_post_wean_TCRB"
assert_columns(repertoire,c(fcb_col))
fig2h_membership <- bind_rows(lapply(names(fig2h_categories),function(cat) {
  ids <- fig2h_categories[[cat]]
  d <- data.frame(Category=cat,CDR3=ids,stringsAsFactors=FALSE)
  left_join(d,repertoire |> transmute(CDR3,
    FcB_productive_frequency=.data[[fcb_col]]),by="CDR3") |>
    mutate(FcB_productive_frequency=coalesce(FcB_productive_frequency,0),
      Detected_in_FcB=FcB_productive_frequency>0,
      # Estimates rounded per clonotype; compare with original immunoSEQ counts.
      Estimated_FcB_templates=round(FcB_productive_frequency*bulk_fcb_total_templates))
}))
fig2h_summary <- fig2h_membership |> group_by(Category) |>
  summarise(Total_clonotypes=n(),Detected_in_FcB=sum(Detected_in_FcB),
    Clonal_detection_rate_percent=100*Detected_in_FcB/Total_clonotypes,
    Cumulative_productive_frequency_in_FcB_percent=100*sum(FcB_productive_frequency),
    Estimated_total_FcB_templates=sum(Estimated_FcB_templates),
    Estimated_templates_per_detected_clonotype=ifelse(Detected_in_FcB>0,
      Estimated_total_FcB_templates/Detected_in_FcB,NA_real_),.groups="drop")
write_table(fig2h_membership,"Fig2H_clone_membership.csv")
write_table(fig2h_summary,"Fig2H_biopsy_enrichment_summary.csv")
# Explicit audit against independently checked five antigen-assigned categories.
expected <- data.frame(Category=c("Pre-DRTCC only","Post-DRTCC only","SRTCC only",
 "Pre-CRTCC","Post-CRTCC"),Total_clonotypes=c(1050,161,248,10,7),
 Detected_in_FcB=c(16,44,59,2,5))
check <- left_join(expected,fig2h_summary,by="Category",suffix=c("_expected","_observed"))
write_table(check,"Audit_Fig2H_expected_vs_observed.csv")
if(any(check$Total_clonotypes_expected!=check$Total_clonotypes_observed) ||
   any(check$Detected_in_FcB_expected!=check$Detected_in_FcB_observed))
 warning("Fig. 2H counts differ from previously validated results; inspect audit CSV")

# --- Section 6: Fig. 2I/J frequency and competition rank ---------------------
# Biopsy ranks are within BULK liver tissue, not sorted CD8 cells.
sample_cols <- c(`Blood 0M`="UPITT_1006_UNSTIM_CD8_TCRB",
 `Blood 6M`="UPITT_1006_6mth_CD8_TCRB",
 `Blood 12M`="UPITT_1006_12mth_CD8_TCRB",
 `Blood 30M`="30M_Unstim_CD8_TCRB",`Blood 54M`="54M_Unstim_CD8_TCRB",
 `MLR 0M`="UPITT_1006_CFSElow_CD8_TCRB",
 `MLR 30M`="30M_CFSE_Low_CD8_TCRB",
 `MLR 33.7M`="33-7M_CFSE_low_CD8_TCRB",
 `MLR 54M`="54M_CFSE_Low_CD8_TCRB",
 `Spike 30M`="30M_SCoV2_S_CD8_TCRB",
 `Spike 33.7M`="33-7M_SCoV2_S_CD8_TCRB",
 `Spike 54M`="54M_SCoV2_S_CD8_TCRB",
 `Liver 12M (bulk)`="UPITT_1006_12mth_LiverBx_TCRB",
 `Liver 33.7M (bulk)`=fcb_col)
assert_columns(repertoire,unname(sample_cols))
if (!all(representative_crtccs %in% all_crtcc))
  stop("One or more Fig. 2I/J representative clones are not in the curated CRTCC set")
rank_table <- data.frame(CDR3=repertoire$CDR3)
for(nm in names(sample_cols)) rank_table[[nm]] <- competition_rank(repertoire[[sample_cols[[nm]]]])
freq_table <- repertoire |> filter(CDR3 %in% representative_crtccs) |>
  select(CDR3,all_of(unname(sample_cols)))
names(freq_table)[-1] <- names(sample_cols)
freq_table[-1] <- lapply(freq_table[-1],function(x) 100*x)
rank_display <- rank_table |> filter(CDR3 %in% representative_crtccs)
fig2ij_long <- freq_table |> pivot_longer(-CDR3,names_to="Sample",
  values_to="Productive_frequency_percent") |>
  left_join(rank_display |> pivot_longer(-CDR3,names_to="Sample",
    values_to="Abundance_rank"),by=c("CDR3","Sample")) |>
  mutate(Rank_display=ifelse(is.na(Abundance_rank),"ND",as.character(Abundance_rank)))
write_table(fig2ij_long,"Fig2IJ_representative_CRTCC_frequency_rank.csv")
wb <- openxlsx::createWorkbook()
for (nm in c("Productive_frequency_percent","Abundance_rank","Long_source_data","Notes"))
  openxlsx::addWorksheet(wb,nm)
openxlsx::writeData(wb,"Productive_frequency_percent",freq_table)
openxlsx::writeData(wb,"Abundance_rank",rank_display)
openxlsx::writeData(wb,"Long_source_data",fig2ij_long)
openxlsx::writeData(wb,"Notes",data.frame(Note=c(
  "Rank 1 is the most abundant detected clonotype in the indicated repertoire.",
  "Competition ranking: ties share rank; undetected clonotypes are ND.",
  "Biopsy ranks use bulk liver-tissue repertoire; MLR and spike use sorted CD8 repertoires.")))
openxlsx::saveWorkbook(wb,file.path(tables_dir,"Fig2IJ_representative_CRTCC_frequency_rank.xlsx"),overwrite=TRUE)

# --- Section 7: reproducibility report ---------------------------------------
write_table(data.frame(Check=c("Pre-DRTCC missing CD4","Pre-DRTCC missing CD8",
 "Complete post-CD8 DRTCC","Curated CD8 SRTCC","Pre-CRTCC","Post-CRTCC",
 "Total CRTCC","Fig2F selected","Comparator selected","Fig2H category overlaps"),
 Value=c(length(missing_pre_cd4),length(missing_pre_cd8),length(post_drtcc_complete),
 length(all_srtcc),length(pre_crtcc),length(post_crtcc),length(all_crtcc),
 nrow(fig2f_wide),length(comparator),sum(overlap_audit$Overlap_count))),
 "Audit_reproducibility_summary.csv")
writeLines(capture.output(sessionInfo()),file.path(results_dir,"sessionInfo.txt"))
message("Finished. Review results/tables/Audit_* before replacing the published analysis.")
