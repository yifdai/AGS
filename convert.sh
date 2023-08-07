#!/bin/bash
VCF_NAME="$1"
VCF_PATH="$2"

dir="/mount/result"

if [ ! -d "${dir}" ]; then 
     mkdir "${dir}"
fi

if [ ! -d "${dir}/${VCF_NAME}" ]; then 
     mkdir "${dir}/${VCF_NAME}"
fi

DIR="$dir/${VCF_NAME}"

echo "Running plink commands"

/opt/plink/plink --vcf "${VCF_PATH}" --recode --make-bed --out "${DIR}/${VCF_NAME}"
/opt/plink/plink --bfile "${DIR}/${VCF_NAME}" --recodeA --out "${DIR}/${VCF_NAME}"



