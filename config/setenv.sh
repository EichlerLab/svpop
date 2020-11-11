# Source this file to load the environment

CENTOS_RELEASE=$(lsb_release -sir | sed -re 's/CentOS\s+([0-9]+).*/\1/')

if [ -z ${CENTOS_RELEASE} ]; then
    CENTOS_RELEASE=Unknown
fi

if [[ ${CENTOS_RELEASE} == 6 ]]; then
    module purge
    module load modules modules-init modules-gs/prod modules-eichler/prod

    module load miniconda/4.5.12
    module load bedtools/2.28.0
    module load htslib/1.9
    module load samtools/1.9
    module load bcftools/1.9
    module load vcflib/20170824
    module load ucsc/20160823
    #module load snpeff/4.3p
    #module load R/3.4.0
    module load perl/5.14.2 ensembl-vep/90.6 RepeatMasker/3.3.0
    module load minimap2/2.17
    module load pbconda/201910
elif [[ ${CENTOS_RELEASE} == 7 ]]; then
    module purge
    module load modules modules-init modules-gs/prod modules-eichler/prod

    module load miniconda/4.5.12
    module load bedtools/2.29.0
    module load htslib/1.9-20
    module load samtools/1.9
    module load bcftools/1.9
    module load vcflib/202002
    module load ucsc/202003
    module load RepeatMasker/4.1.0
    module load minimap2/2.17
    module load pbconda/201910
else
    echo "Warning: Unknown CentOS release: ${CENTOS_RELEASE}: No modules loaded"
fi

set -euo pipefail
