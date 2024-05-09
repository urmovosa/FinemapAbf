

set -f

nextflow_path=[folder with Nextflow executable]

inputDir=[Folder containing GWAS sumstats in .parquet format]
outputDir=[Output folder]
referenceDir=[Folder containing SNP reference file]

NXF_VER=21.10.6 ${nextflow_path}/nextflow run main.nf  \
--inputDir ${inputDir} \
--outputDir ${outputDir} \
--SnpRefFile ${referenceDir} \
-resume \
-profile slurm,singularity
