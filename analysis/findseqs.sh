### A script to select which sequecnes mataches our primers

## Usage bash findseqs.sh <folder with seqs> <primerF> <primerR>

## Before launch activate your cutadapt environment

folder=${1}

output_folder="${folder}"/trimmed

FWD=${2}

REV=${3}

for target in "${folder}"/*.fasta ; do
base_target="$(basename ${target})"
name="${base_target%.*}"

cutadapt -g ""${FWD}";min_overlap=15" --discard-untrimmed -e 0.2 "${target}" --rc \
-o "${output_folder}"/"${name}"_with_fwd.fasta > "${output_folder}"/"${name}".report_1

cutadapt -g ""${REV}";min_overlap=15" --discard-untrimmed -e 0.2 "${output_folder}"/"${name}"_with_fwd.fasta --rc \
-o "${output_folder}"/"${name}"_with_both.fasta > "${output_folder}"/"${name}".report_2

done