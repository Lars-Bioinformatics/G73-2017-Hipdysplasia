# install.packages("vcfR")
library(vcfR)
library(tidyverse)

# setwd("~/OneDrive - Syddansk Universitet/Projekter/G73-2017-Hoftedysplasi/")
setwd("/work/sduvarcall/G73-2017-Hoftedysplasi/")

# files = list.files("vcf_gnomAD/vcf_gnomAD_MAF_filtered/")
files = list.files("vcf_gnomAD/vcf_gnomAD_annotated/")
file = files[1]

# Output
dir.create("VAF_distribution", showWarnings = F)

# Plot VAF density
for (file in files){
  # vcf = read.vcfR("vcf_gnomAD/vcf_gnomAD_MAF_filtered/fam1_variants_gnomAD_filterMAF_onePercent.vcf.gz")
  # vcf = read.vcfR(paste0("vcf_gnomAD/vcf_gnomAD_MAF_filtered/",file))
  vcf = read.vcfR(paste0("vcf_gnomAD/vcf_gnomAD_annotated/",file))
  fam = str_split(file, "_", n = 2, simplify = T)[1]
  
  # vcf@gt
  df = vcfR::extract_gt_tidy(vcf, format_fields = c("AD","DP", "GT"))
  # df = vcfR::extract_gt_tidy(vcf)
  readCounts = df %>%
    # Keep only "simple" genotypes i.e. 0/1, 1/1
    filter(gt_GT %in% c("0/1","0|1","1/1","1|1")) %>%
    # Only allow variants with one type of variant
    filter(str_count(gt_AD,",")==1) %>%
    # Only keep SNVs
    filter(str_length(gt_GT_alleles)==3) %>%
    # Split RefDepth and AlleleDepth
    separate(gt_AD, c("RD","AD"), sep = ",", convert = T) %>%
    # filter(AD>0) %>%
    # Compute VAF
    mutate(VAF=AD/gt_DP)
  
  p = readCounts %>%
    ggplot(aes(x=VAF, color=Indiv, fill=Indiv)) +
    geom_density(alpha=0.5) +
    theme_minimal() +
    ggtitle(paste("Family", str_remove(fam,"fam"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~Indiv)
  p
  ggsave(filename = paste0("VAF_distribution/", fam, "_VAF_distribution.png"),
         plot = p, width=12, height = 8, bg = "white")
  
  
}

# Het estimates
het_estimates = map_dfr(
  .x = files,
  .f = function(file){
    vcf = read.vcfR(paste0("vcf_gnomAD/vcf_gnomAD_annotated/",file))
    fam = str_split(file, "_", n = 2, simplify = T)[1]
    
    # vcf@gt
    df = vcfR::extract_gt_tidy(vcf, format_fields = c("AD","DP", "GT"))
    # df = vcfR::extract_gt_tidy(vcf)
    readCounts = df %>%
      # Keep only "simple" genotypes i.e. 0/1, 1/1
      filter(gt_GT %in% c("0/1","0|1","1/1","1|1")) %>%
      # Only allow variants with one type of variant
      filter(str_count(gt_AD,",")==1) %>%
      # Only keep SNVs
      filter(str_length(gt_GT_alleles)==3)
    
    # Het estimate
    het_estimate = readCounts %>% 
      mutate(gt=ifelse(gt_GT %in% c("0/1", "0|1"),"het_count","hom_count")) %>% 
      count(Indiv, gt) %>%
      pivot_wider(names_from = gt, values_from = n) %>%
      mutate(het_ratio=het_count/(het_count+hom_count), family=fam) %>%
      select(sample=Indiv, family, het_count, hom_count, het_ratio)
    return(het_estimate)
  }
)
het_estimates
dir.create("Heterozygous_ratios/", showWarnings = F)
write.table(het_estimates, file = "Heterozygous_ratios/Heterozygous_ratios.txt", quote = F, row.names = F, sep = "\t")
