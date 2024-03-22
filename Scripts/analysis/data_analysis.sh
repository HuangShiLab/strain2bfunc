feature_table=$1
metadata=$2
out_dir=$3
dist_file=$4
prefix=$5

# Check if the file exists
if [ ! -f "$feature_table" ]; then
    echo "File $feature_table does not exist."
    exit 1
fi

if [ ! -f "$metadata" ]; then
    echo "File $metadata does not exist."
    exit 1
fi

# Continue with the rest of the script

#Abundance distribution visualization
mkdir -p ${out_dir}/Abundance_distribution
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Distribution.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Abundance_distribution -p ${prefix}"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Distribution.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Abundance_distribution -p ${prefix}

#Multivariate Statistical Analysis
mkdir -p ${out_dir}/Alpha_Diversity
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Adiversity.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Alpha_Diversity -p ${prefix}"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Adiversity.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Alpha_Diversity -p ${prefix}

#Dist Calculation
mkdir -p ${out_dir}/Distance
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Dist.R -i ${feature_table} -o ${out_dir}/Distance/${dist_file}"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Dist.R -i ${feature_table} -o ${out_dir}/Distance/${dist_file}
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Heatmap.R -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Distance/dist_heatmap.pdf"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Heatmap.R -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Distance/dist_heatmap.pdf
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Hcluster.R -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Distance/dist_cluster.pdf"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Hcluster.R -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Distance/dist_cluster.pdf

# Get the number of lines in the feature table
line_count=$(wc -l < ${feature_table})

# Check if the number of lines is less than or equal to 4
mkdir -p ${out_dir}/Clustering/
if [ "$line_count" -le 4 ]; then
    echo "The number of strains of the samples under this species is not greater than 3, Unable to perform PCoA and PCA analysis." 
    echo "The number of strains of the samples under this species is not greater than 3, Unable to perform PCoA and PCA analysis." > ${out_dir}/Clustering/Clustering_log.txt
else
    #PCoA Calculation
    echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Pcoa.R -m ${metadata} -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Clustering/PCoA.pdf"
    Rscript ${Strain2bFunc}/Scripts/analysis/PM_Pcoa.R -m ${metadata} -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Clustering/PCoA.pdf

    #PCA Calculation
    echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Pca.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Clustering/PCA.pdf"
    Rscript ${Strain2bFunc}/Scripts/analysis/PM_Pca.R -m ${metadata} -i ${feature_table} -o ${out_dir}/Clustering/PCA.pdf
fi

mkdir -p ${out_dir}/Beta_Diversity
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Bdiversity.R -m ${metadata} -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Beta_Diversity -p ${prefix} -n \"\" -t 8"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Bdiversity.R -m ${metadata} -d ${out_dir}/Distance/${dist_file} -o ${out_dir}/Beta_Diversity -p ${prefix} -n "" -t 8

#Marker Analysis
mkdir -p ${out_dir}/Markers
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_Test.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix} -P F"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_Test.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix} -P F
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_RFscore.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix}"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_RFscore.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix}
echo "Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_Corr.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix}"
Rscript ${Strain2bFunc}/Scripts/analysis/PM_Marker_Corr.R -i ${feature_table} -m ${metadata} -o ${out_dir}/Markers -p ${prefix}

