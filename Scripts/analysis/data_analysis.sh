feature_table=$1
metadata=$2
out_dir=$3
dist_file=$4
prefix=$5

#Dist Calculation
mkdir -p ${out_dir}/Distance
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Dist.R -i ${feature_table} -o ${out_dir}/Distance/${dist_file}

#PCoA Calculation
mkdir -p ${out_dir}/Clustering/
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Pcoa.R -m ${metadata} -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Clustering/PCoA.pdf

#PCA Calculation
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Pca.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Clustering/PCA.pdf

#Multivariate Statistical Analysis
mkdir -p ${out_dir}/Alpha_Diversity
mkdir -p ${out_dir}/Beta_Diversity
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Adiversity.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Alpha_Diversity -p ${prefix}
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Bdiversity.R -m ${metadata} -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Beta_Diversity -p ${prefix} -n "" -t 8

#Marker Analysis
mkdir -p ${out_dir}/Markers
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_Test.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix} -P F
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_RFscore.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix}
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_Corr.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix}

