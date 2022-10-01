for r in {"R1","R2"}; do
    echo "cat G73-2874_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz G73-2874_nimblegen-medexome_HM5TKBGX5_${r}.fastq.gz > G73-2874_nimblegen-medexome_${r}.fastq.gz"
    cat G73-2874_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz G73-2874_nimblegen-medexome_HM5TKBGX5_${r}.fastq.gz > G73-2874_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-2874_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz G73-2874_nimblegen-medexome_HM5TKBGX5_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-2A15_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-2A15_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz > G73-2A15_nimblegen-medexome_${r}.fastq.gz"
    cat G73-2A15_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-2A15_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz > G73-2A15_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-2A15_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-2A15_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-5504_nimblegen-medexome_HC3C5BGX9_${r}.fastq.gz G73-5504_nimblegen-medexome_HTJ2YBGX5_${r}.fastq.gz > G73-5504_nimblegen-medexome_${r}.fastq.gz"
    cat G73-5504_nimblegen-medexome_HC3C5BGX9_${r}.fastq.gz G73-5504_nimblegen-medexome_HTJ2YBGX5_${r}.fastq.gz > G73-5504_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-5504_nimblegen-medexome_HC3C5BGX9_${r}.fastq.gz G73-5504_nimblegen-medexome_HTJ2YBGX5_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-5C3C_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-5C3C_nimblegen-medexome_HMHV5BGX7_${r}.fastq.gz > G73-5C3C_nimblegen-medexome_${r}.fastq.gz"
    cat G73-5C3C_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-5C3C_nimblegen-medexome_HMHV5BGX7_${r}.fastq.gz > G73-5C3C_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-5C3C_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-5C3C_nimblegen-medexome_HMHV5BGX7_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-96D6_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-96D6_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz > G73-96D6_nimblegen-medexome_${r}.fastq.gz"
    cat G73-96D6_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-96D6_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz > G73-96D6_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-96D6_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-96D6_nimblegen-medexome_HKGY5BGX7_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-CE38_nimblegen-medexome_HM5TKBGX5_${r}.fastq.gz G73-CE38_nimblegen-medexome_HYFYTBGX5_${r}.fastq.gz > G73-CE38_nimblegen-medexome_${r}.fastq.gz"
    cat G73-CE38_nimblegen-medexome_HM5TKBGX5_${r}.fastq.gz G73-CE38_nimblegen-medexome_HYFYTBGX5_${r}.fastq.gz > G73-CE38_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-CE38_nimblegen-medexome_HM5TKBGX5_${r}.fastq.gz G73-CE38_nimblegen-medexome_HYFYTBGX5_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-E6D8_nimblegen-medexome_HC3C5BGX9_${r}.fastq.gz G73-E6D8_nimblegen-medexome_HMHV5BGX7_${r}.fastq.gz > G73-E6D8_nimblegen-medexome_${r}.fastq.gz"
    cat G73-E6D8_nimblegen-medexome_HC3C5BGX9_${r}.fastq.gz G73-E6D8_nimblegen-medexome_HMHV5BGX7_${r}.fastq.gz > G73-E6D8_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-E6D8_nimblegen-medexome_HC3C5BGX9_${r}.fastq.gz G73-E6D8_nimblegen-medexome_HMHV5BGX7_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-E8BE_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-E8BE_nimblegen-medexome_HTJ2YBGX5_${r}.fastq.gz > G73-E8BE_nimblegen-medexome_${r}.fastq.gz"
    cat G73-E8BE_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-E8BE_nimblegen-medexome_HTJ2YBGX5_${r}.fastq.gz > G73-E8BE_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-E8BE_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-E8BE_nimblegen-medexome_HTJ2YBGX5_${r}.fastq.gz multi_fastq_samples
    echo "cat G73-FDD5_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-FDD5_nimblegen-medexome_HVWL7BGX5_${r}.fastq.gz > G73-FDD5_nimblegen-medexome_${r}.fastq.gz"
    cat G73-FDD5_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-FDD5_nimblegen-medexome_HVWL7BGX5_${r}.fastq.gz > G73-FDD5_nimblegen-medexome_merged_${r}.fastq.gz
    mv G73-FDD5_nimblegen-medexome_HJ55FBGX9_${r}.fastq.gz G73-FDD5_nimblegen-medexome_HVWL7BGX5_${r}.fastq.gz multi_fastq_samples
done