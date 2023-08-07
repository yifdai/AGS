#!/bin/bash
VCF_NAME="$1"
VCF_PATH="$2"

echo "$VCF_PATH"
echo "$VCF_NAME"

dir="/mount/result"

if [ ! -d "${dir}/${VCF_NAME}/rrblup-result" ]; then
  mkdir "${dir}/${VCF_NAME}/rrblup-result"
fi

DIR="${dir}/${VCF_NAME}"

/opt/plink/plink --bfile "${DIR}/${VCF_NAME}" --make-grm-bin  --allow-extra-chr --out "${DIR}/${VCF_NAME}"

