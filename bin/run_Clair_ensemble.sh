#!/usr/bin/env bash
while getopts b:d:r:c:m:e:t:o: option
do
    case "${option}"
        in
        # assumed given bam file is sorted and indexed
        b) BAM_FILE_PATH=`readlink -f ${OPTARG}`;;
        d) BED_FILE_PATH=`readlink -f ${OPTARG}`;;
        r) REFERENCE_FASTA_FILE_PATH=`readlink -f ${OPTARG}`;;
        c) CLAIR=`readlink -f ${OPTARG}`;;
        m) old_IFS=$IFS
           IFS=','
           CLAIR_MODELS=($OPTARG)
           IFS=${old_IFS};;
        e) ENSEMBLE_CPP_EXECUTABLE=`readlink -f ${OPTARG}`;;
        t) PARALLEL_THREADS=${OPTARG};;
        o) ROOT_FOLDER_PATH=`readlink -f ${OPTARG}`;;
    esac
done

SCRIPT_DIR_PATH=`dirname readlink -f $0`
NO_OF_CLAIR_MODELS=${#CLAIR_MODELS[@]}
DATE_TIME=`date "+%Y%m%d_%H%M%S"`
WORKING_DIRECTORY="${ROOT_FOLDER_PATH}/${DATE_TIME}"
mkdir -p ${WORKING_DIRECTORY}
BAM_FILE_PATHS=("${BAM_FILE_PATH}")
INTERMEDIATE_OUTPUT_FOLDER="${WORKING_DIRECTORY}/tmp_output"
mkdir  ${INTERMEDIATE_OUTPUT_FOLDER}
cd ${INTERMEDIATE_OUTPUT_FOLDER}
for i in "${!BAM_FILE_PATHS[@]}"
do
  INPUT_BAM_FILE_PATH="${BAM_FILE_PATHS[i]}"
  SAMPLE_NAME="MES"
  BAM_PREFIX=`printf "%02d" ${i}`

  for j in "${!CLAIR_MODELS[@]}"
  do
    CLAIR_MODEL="${CLAIR_MODELS[j]}"
    MODEL_PREFIX=`printf "%02d" ${j}`
    SCRIPT_OUTPUT_FOLDER="m${MODEL_PREFIX}_b${BAM_PREFIX}"

    mkdir ${SCRIPT_OUTPUT_FOLDER}
    OUTPUT_PREFIX="${SCRIPT_OUTPUT_FOLDER}/tmp"

    python ${CLAIR} callVarBamParallel \
    --chkpnt_fn "${CLAIR_MODEL}" \
    --ref_fn "${REFERENCE_FASTA_FILE_PATH}" \
    --bed_fn "${BED_FILE_PATH}" \
    --bam_fn "${INPUT_BAM_FILE_PATH}" \
    --pysam_for_all_indel_bases \
    --output_for_ensemble \
    --refChunkSize "10000000" \
    --sampleName "${SAMPLE_NAME}" \
    --output_prefix "${OUTPUT_PREFIX}" > ${SCRIPT_OUTPUT_FOLDER}/call.sh
  done
done
( time cat */call.sh | parallel -j${PARALLEL_THREADS} ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.call.txt

FILES=(`ls m00_b00/*.vcf`)
ENSEMBLE_OUTPUT_FOLDER="${INTERMEDIATE_OUTPUT_FOLDER}/ensemble"
mkdir ${ENSEMBLE_OUTPUT_FOLDER}
MININUM_NO_OF_VOTE_FOR_VARIANT="$(((${#BAM_FILE_PATHS[@]}*${#CLAIR_MODELS[@]}+2)/2))"
rm ensemble_command.sh
for i in "${!FILES[@]}"
do
  TARGET_FILE_NAME="${FILES[i]:8}"
  CAT_COMMAND=""
  for j in "${!BAM_FILE_PATHS[@]}"
  do
    BAM_PREFIX=`printf "%02d" ${j}`
    for k in "${!CLAIR_MODELS[@]}"
    do
      MODEL_PREFIX=`printf "%02d" ${k}`
      FOLDER_NAME="m${MODEL_PREFIX}_b${BAM_PREFIX}"
      CAT_COMMAND="${CAT_COMMAND} ${FOLDER_NAME}/${TARGET_FILE_NAME}"
    done
  done


  echo "cat ${CAT_COMMAND:1} | ${ENSEMBLE_CPP_EXECUTABLE} ${MININUM_NO_OF_VOTE_FOR_VARIANT} > ${ENSEMBLE_OUTPUT_FOLDER}/${TARGET_FILE_NAME}" >> ensemble_command.sh
done
( time cat ensemble_command.sh | parallel -j${PARALLEL_THREADS} ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.ensemble.txt


VCF_OUTPUT_FOLDER="${WORKING_DIRECTORY}/output"
mkdir ${VCF_OUTPUT_FOLDER}
cd ${WORKING_DIRECTORY}
INPUT_FILES=(`ls tmp_output/ensemble/*.vcf`)

for i in "${!INPUT_FILES[@]}"
do
  FILE_NAME="${INPUT_FILES[i]:20}"
  echo "cat tmp_output/ensemble/${FILE_NAME} | \
  python ${CLAIR} call_var \
  --chkpnt_fn "${CLAIR_MODEL}" \
  --ref_fn "${REFERENCE_FASTA_FILE_PATH}" \
  --bam_fn "${BAM_FILE_PATH}" \
  --call_fn "${VCF_OUTPUT_FOLDER}/${FILE_NAME}" \
  --sampleName "MES" \
  --pysam_for_all_indel_bases \
  --input_probabilities" >> output.sh

done
( time cat output.sh | parallel -j${PARALLEL_THREADS} ) |& tee ${INTERMEDIATE_OUTPUT_FOLDER}/log.output.txt


cd ${WORKING_DIRECTORY}
vcfcat ${VCF_OUTPUT_FOLDER}/*.vcf | sort -k1,1V -k2,2n > snp_and_indel.vcf
cat snp_and_indel.vcf | python ${SCRIPT_DIR_PATH}/Clair-ensemble/Clair.beta.ensemble.cpu/clair/post_processing/overlap_variant.py > snp_and_indel.filtered.vcf



