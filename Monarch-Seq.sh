#!/bin/bash

set -e

VERSION=0.55
# Log message function (normal text) with Step
log() {
        echo -ne "[ $(date +%m-%d-%y" "%H:%M:%S) | $1 ] $2"
}
# Warning message function (yellow text)
warn() {
        echo -e "\033[1;33m[ $(date +%m-%d-%y" "%H:%M:%S) | $1 ] $2\033[0m"
        }
# "Done" message function
DONEmsg() {
        echo -e "[ $(date +%m-%d-%y" "%H:%M:%S) | $1 ] Done!"
}
# "Bye" message function (green text)
BYEmsg() {
        echo -e "\033[1;32m[ $(date +%m-%d-%y" "%H:%M:%S) ] ALL DONE!"
        echo -e "
This is a goodbye message.
Probably will have citation, mail contact, and more.
Bye\033[0m"
}
# thanks, stackoverflow!
join_path() {
    echo "${1:+$1/}$2" | sed 's#//#/#g'
}

# usage function
usage() {
echo -e "
Monarch-Seq
Automated pipeline to find circular RNAs using Circ-Seq reads
----------------------------------------------------------------

\e[4mUsage\e[0m:

 Monarch-Seq.sh --ref_genome <reference.fasta> --reads <reads1.fasta,reads2.fasta>
                --output <path/to/output> [--threads INT]
               [--prefix STR] [--max_circle_size INT]
               [--min_sec_hit_size INT] [--min_read_cov FLOAT]
               [--max_mismatch INT] [--no_strandness] [--dont_invert_strand]

\e[4mRequired arguments\e[0m:

  --ref_genome         The reference genome to be used as subject for blast

  --reads              Reads to be used as query for blast. Or a comma separated
                       list of reads, if applicable

  --output             Path to output dir. Will be created if does not exist

\e[4mOptional arguments\e[0m:

 --threads             Number of threads for the analysis (Default = 10)

 --prefix              Prefix for files created by Monarch-Seq (Default = circles)

 --max_circle_size     Maximum circle size (Default = 3500)

 --min_sec_hit_size    Minimum size of the second hit (Default = 8)

 --min_read_cov        Minimum read coveraged by first and second hits together.
                       Must be number between 0 and 1 (0 - 100%) (Default = 0.9)

 --max_mismatch        Maximum number of mismatches allowed for each alignment
                       (Default = 0)

 --no_strandness        Do not consider sequencing protocol stranded. Monarch-Seq
                       default behavior is to consider sequencing protocol stranded.
                       This paramenter overrides --dont_invert_strand.

 --dont_invert_strand  Do not invert input strand. Default behavior of Monarch-Seq
                       is to turn 'plus' strand turns into 'minus' and vice versa.
                       Depending on the sequencing protocol one might not want this.
"

exit 2
}

# error function
error_exit() {
        msg=$1
        exit_code=${2:-1}
        echo -e "\033[1;31m${msg}\033[0m"
        usage
}

check_deps(){
	for app in makeblastdb blastn parallel; do
                        command -v $app >/dev/null 2>&1 || \
                        error_exit "ERROR: Cannot find ${app} in your PATH variable"
        done
}

# option strings
SHORT=h
LONG=help,ref_genome:,reads:,prefix:,max_circle_size:,output:,min_sec_hit_size:,threads:,min_read_cov:,max_mismatch:,dont_invert_strand,no_strandness

# read the options
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")
if [ $? != 0 ] ; then echo "Failed to parse options...exiting." ; exit 1 ; fi
eval set -- "$OPTS"

# set arguments initial values
THREADS=10
DELTA=3
MAX_CIRCLE_SIZE=3500
MIN_READ_COV=0.9
PREFIX="circles"
MAX_MISMATCH=0
MIN_SEC_HIT_SIZE=8
INVERT_STRAND="YES"
STRANDNESS="YES"

# extract options and their arguments into variables.
while true ; do
  case "$1" in
    -h | --help )
      usage
      ;;
    --ref_genome )
      REF_GENOME="$2"
      shift 2
      ;;
    --threads )
      THREADS="$2"
      shift 2
      ;;
    --reads )
      READS="$2"
      shift 2
      ;;
    --max_circle_size )
      MAX_CIRCLE_SIZE="$2"
      shift 2
      ;;
    --min_sec_hit_size )
      MIN_SEC_HIT_SIZE="$2"
      shift 2
      ;;
    --output )
      OUTPUT="$2"
      shift 2
      ;;
    --min_read_cov )
      MIN_READ_COV="$2"
      shift 2
      ;;
    --prefix )
      PREFIX="$2"
      shift 2
      ;;
    --max_mismatch )
      MAX_MISMATCH="$2"
      shift 2
      ;;
    --dont_invert_strand )
      INVERT_STRAND="NO"
      shift
      ;;
    --no_strandness )
      STRANDNESS="NO"
      shift
        ;;
    -- )
      shift
      break
      ;;
    *)
      error_exit "ERROR: Please, supply required arguments correctly."
      ;;
  esac
done

# setting home dir, for utils calling
tmp=$(realpath "$0")
HOME_DIR=${tmp%/*}

check_deps

# checking if any reuquired args is empty
if [ -z "${REF_GENOME}" ] || [ -z "${READS}" ] || [ -z "${OUTPUT}" ]; then
        error_exit "ERROR: Please, supply required arguments correctly."
fi

# Starting Monarch-Seq
echo "=========================================================================
                           Monarch-Seq (v$VERSION)
      Automated pipeline to find circular RNAs using Circ-Seq reads
-------------------------------------------------------------------------
REFERENCE GENOME: $REF_GENOME
READS: $READS
OUTPUT: $OUTPUT
PREFIX: $PREFIX
THREADS: $THREADS
MAX CIRCLE SIZE: $MAX_CIRCLE_SIZE
MIN READ COVERAGE: $MIN_READ_COV
MAX MISMATCH ALLOWED: $MAX_MISMATCH
MIN SECOND HIT SIZE: $MIN_SEC_HIT_SIZE
STRANDED PROTOCOL: $STRANDNESS
INVERT STRAND? $INVERT_STRAND
=========================================================================
"

mkdir -p "$OUTPUT"
if [[ ! -f "$OUTPUT"/reads00.fa && ! -f "$OUTPUT"/"$PREFIX".blastn.tbl.gz ]]; then
        # ================
        # ==== STEP 1 ====
        # ================

        log "Step1" "Creating BLAST database\n"
        makeblastdb -in "$REF_GENOME" -dbtype nucl > /dev/null
        DONEmsg "Step1"

        # ================
        # ==== STEP 2 ====
        # ================
        # concatenate reads, if applicable
        # 1) Cat treated files
        cat_treated=$(mktemp)
        #cat_treated=${OUTPUT}/cat_treated ; touch $cat_treated
        for file in $(echo $READS | sed 's/,/ /g'); do
        	cat $file >> $cat_treated
        done

        splitted=$(join_path "$OUTPUT" reads)
	log "Step2" "Splitting input reads\n"
	split -l 10000 -d --additional-suffix=".fa" "$cat_treated" "$splitted"
	DONEmsg "Step2"
else
        warn "Step1" "Blast database already exist. Skipping..."
        warn "Step2" "Splitted files already exist. Skipping..."
fi


# ================
# ==== STEP 3 ====
# ================
N_FILES=$(find "$OUTPUT" -name "*.fa" | wc -l)
if [[ ! -f "$OUTPUT"/reads00.fa.out && ! -f "$OUTPUT"/"$PREFIX".blastn.tbl.gz ]]; then
	log "Step3" "BLASTn\n"
	AUX=0
	for input in $(find "$OUTPUT" -name "*.fa"); do
		file=$(basename $input)
		sem --will-cite --id $$ --max-procs $THREADS \
			blastn -num_threads 1 -db $REF_GENOME -query $input \
			-outfmt \"6 sseqid qseqid qlen length qstart qend sstart send sstrand mismatch\" \
			-out "$input".out -ungapped -task  blastn-short
		AUX=$(( AUX + 1 ))
		log "Step3" "Running part $AUX of $N_FILES...\r"
	done
        sem --will-cite --id $$ --wait
	echo -ne "\n"
	DONEmsg "Step3"
else
	warn "Step3" "BLASTn already run. Skipping..."
fi

# ================
# ==== STEP 4 ====
# ================

if [[ ! -f "$OUTPUT"/reads00.fa.out.bed && ! -f "$OUTPUT"/"$PREFIX".blastn.tbl.gz ]]; then
        log "Step4" "Parsing Results\n"
        AUX=0
        for input in $(find "$OUTPUT" -name "*.fa.out"); do
		sem --will-cite --id $$ --max-procs $THREADS \
			"$CONDA_PREFIX"/bin/python3 \
			"$HOME_DIR"/scripts/Monarch-GetJunctions.py "$input" $DELTA \
			$MAX_CIRCLE_SIZE $MIN_SEC_HIT_SIZE $MIN_READ_COV \
			$MAX_MISMATCH $INVERT_STRAND $STRANDNESS
		AUX=$(( AUX + 1 ))
		log "Step4" "Running part $AUX of $N_FILES...\r"
	done
	sem --will-cite --id $$ --wait
	echo -ne "\n"
	DONEmsg "Step4"
else
	warn "Step4" "BLAST results already parsed. Skipping..."
fi

# ================
# ==== STEP 5 ====
# ================

if [[ ! -f "$OUTPUT"/"$PREFIX".junctions.info && ! -f "$OUTPUT"/"$PREFIX".blastn.tbl.gz ]]; then
	log "Step5" "Merging files\n"
	#cat ${OUTPUT}/*.fa.out | gzip > "$OUTPUT"/"$PREFIX".blastn.tbl.gz
 	find "$OUTPUT" -name "*.fa.out" -exec cat '{}' \+ | pigz -p $THREADS > "$OUTPUT"/"$PREFIX".blastn.tbl.gz
	#cat ${OUTPUT}/*out.bed > "$OUTPUT"/"$PREFIX".junctions.bed
 	find "$OUTPUT" -name "*out.bed" -exec cat '{}' \+ > "$OUTPUT"/"$PREFIX".junctions.bed
        echo -e "#ReadName\tStrand\tHalf\tSecHitCount\tOverlap\
\tReadStart1\tReadEnd1*\tReadStart2*\tReadEnd2\
\tGenomeStart1\tGenomeEnd1*\tGenomeStart2*\tGenomeEnd2" \
        > "$OUTPUT"/"$PREFIX".junctions.info
        #cat ${OUTPUT}/*.fa.out.info >> "$OUTPUT"/"$PREFIX".junctions.info
	find "$OUTPUT" -name "*.fa.out.info" -exec cat '{}' \+ >> "$OUTPUT"/"$PREFIX".junctions.info
	DONEmsg "Step5"
else
	warn "Step5" "Output files already merged. Skipping..."
fi

# ================
# ==== STEP 6 ====
# ================

# Functions to be used in post processing (To add transcript and junction seqs)
# They are almost the same, but the first is to be used when INVERT_STRAND == YES
#    and the second when INVERT_STRAND == NO.


if [ ! -f "$OUTPUT"/"$PREFIX".junctions.sorted_for_ensembles ]; then
	log "Step6" "Post processing\n"
        log "Step6" "Getting junction sequences\n"
        INFOTMP=$(mktemp)
        AUX=0
        if [ "$INVERT_STRAND" == "YES" ]; then
                # replace 'minus' with '+' and 'plus' with '-'
                cat "$OUTPUT"/"$PREFIX".junctions.info | \
                        sed 's/minus/+/g;s/plus/-/g' | \
                        tail -n +2 | sed 's/\t/franceschini-santos/g' > $INFOTMP
                N_FILES=$(wc -l < $INFOTMP)
                while read -r LINE; do
                        $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs \
                        "$THREADS" ${HOME_DIR}/scripts/Monarch-GetJunctionSeqs.sh \
                                "$LINE" "$OUTPUT" "$REF_GENOME" "$cat_treated" \
                                "$MIN_SEC_HIT_SIZE" "$INVERT_STRAND"
                        AUX=$(( AUX + 1 ))
		        log "Step6" "Running junction $AUX of $N_FILES...\r"
                done < $INFOTMP
                sem --will-cite --id $$ --wait
	        echo -ne "\n"
        else
                # replace 'minus' with '-' and 'plus' with '+'
                cat "$OUTPUT"/"$PREFIX".junctions.info | \
                        sed 's/plus/+/g;s/minus/-/g' | \
                        tail -n +2 | sed 's/\t/franceschini-santos/g' > $INFOTMP
                N_FILES=$(wc -l < $INFOTMP)
                while read -r LINE; do
                        $CONDA_PREFIX/bin/sem --will-cite --id $$ --max-procs \
                        "$THREADS" ${HOME_DIR}/scripts/Monarch-GetJunctionSeqs.sh \
                                "$LINE" "$OUTPUT" "$REF_GENOME" "$cat_treated" \
                                "$MIN_SEC_HIT_SIZE" "$INVERT_STRAND"
                        AUX=$(( AUX + 1 ))
		        log "Step6" "Running junction $AUX of $N_FILES...\r"
                done < $INFOTMP
                sem --will-cite --id $$ --wait
	        echo -ne "\n"
        fi
else
	warn "Step6" "Junction sequences already gotten. Skipping..."
fi
if [ ! -f "$OUTPUT"/"$PREFIX".ensembles.full.readseqs ]; then

    # sorting by chr, start and end
	cat "$OUTPUT"/"$PREFIX".junctions.bed | sort -k1n -k2n -k3n > "$OUTPUT"/"$PREFIX".junctions.sorted_for_ensembles

	# running PostProcessing.py
        log "Step6" "Running ensembles annotation\n"
	"$CONDA_PREFIX"/bin/python3 "$HOME_DIR"/scripts/Monarch-PostProcessing.py \
		"$OUTPUT"/"$PREFIX".junctions.sorted_for_ensembles "$OUTPUT" \
		"$PREFIX" "$cat_treated"
	# writing temp ensemble ID to files .full
	#cat "$OUTPUT"/"$PREFIX".ensembles.full.readnames.tmp \
		#| awk -v OFS="\t" '{ print $1,$2,$3,$4,$1"__"$5,$6,$7,$8,$9,$10 }' \
		#> "$OUTPUT"/"$PREFIX".ensembles.full.readnames
        echo -ne "#CircEnse\tChr\tStart\tEnd\tReadID\tFlag\tStrand\tReadSeq\tTranscriptSize\
\tHalfA\tHalfB\tEmpiricalJunction\tReadlJunction\tTranscript\n" > "$OUTPUT"/"$PREFIX".ensembles.full
	cat "$OUTPUT"/"$PREFIX".ensembles.full.readseqs.tmp \
		| awk -v OFS="\t" '{$5=$1"__"$5; print $0 }' \
		>> "$OUTPUT"/"$PREFIX".ensembles.full
        # Adding transcript and junction sequences to junctions.seqs file

        echo ""
        DONEmsg "Step6"
else
	warn "Step6" "Post processing already done. Skipping..."
fi

# ================
# ==== STEP 7 ====
# ================

if [ -f "$OUTPUT"/reads00.fa.out.bed ]; then
	log "Step7" "Removing temporary files\n"
 	rm -rf ${OUTPUT}/tmp/
     	rm -f "$OUTPUT"/"$PREFIX".junctions.sorted_for_ensembles
	rm -f "$OUTPUT"/"$PREFIX".ensembles.full.read*.tmp
     	find "${OUTPUT}" -maxdepth 1 -type f -name "*.fa" -delete
     	find "${OUTPUT}" -maxdepth 1 -type f -name "*.fa.out.bed" -delete
        find "${OUTPUT}" -maxdepth 1 -type f -name "*.fa.out.info" -delete
	find "${OUTPUT}" -maxdepth 1 -type f -name "*.fa.out" -delete
	DONEmsg "Step7"
else
	warn "Step7" "Temporary files already removed. Skipping..."
fi

BYEmsg
