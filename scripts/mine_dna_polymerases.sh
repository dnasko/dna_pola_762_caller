#!/bin/bash
set -e

usage() {
    echo; echo "Usage: $0 --query=/Path/to/infile.fasta --database=/Path/to/PolA_DB --out=/Path/to/output.fasta [--threads=4]"
    echo "  --fasta     Path to the input AA FASTA file"
    echo "  --out       Path to the output AA FASTA file"
    echo "  --database  Path to the PolA database FASTA file"
    echo "  --threads   Number of threads to use (Default=1)"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

## ========================
## Globals
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
THREADS=1

## ========================
## Process the input parameters
if [ $# -gt 5 ] || [ $# -lt 4 ]
then
    usage
fi

while true
do
    case $1 in
    --help|-h)
	    usage
	    exit;;
    --fasta=?*)
	    FASTA=${1#*=};;
    --fasta|fasta=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --out=?*)
	    OUT=${1#*=};;
    --out|out=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --threads=?*)
	    THREADS=${1#*=};;
    --threads|threads=)
	    echo "$0: missing argument for '$1' option"
	    usage
	    exit 1;;
    --database=?*)
            DATABASE=${1#*=};;
    --database|database=)
            echo "$0: missing argument for '$1' option"
            usage
            exit 1;;
    --)
	    shift
	    break;;
    -?*)
	    echo "$0: invalid option: $1"
	    usage
	    exit 1;;
    *)
	    break
    esac
    shift
done

##===
## Logging
LOG="$( cd "$( dirname "${OUT}" )" && pwd )"
LOG="${LOG}/mine_polas.log"

echo -n "# Building a BLASTp database of the query file ......" | tee -a ${LOG}
if [ -f "${FASTA}.phr" ]; then

else
    makeblastdb -in ${FASTA} -dbtype prot
    status=$?
fi
echo -n " DONE "  | tee -a ${LOG}; date '+%H:%M:%S %Y-%m-%d' |tee -a ${LOG}

echo -n "# BLASTp search for candidiate PolA sequences ......." | tee -a ${LOG}
blastp -query "${DATABASE}" \
       -db "${FASTA}" \
       -out "${OUT}.btab" \
       -num_threads "${THREADS}" \
       -evalue 1e-5 \
       -max_target_seqs 1000000 \
       -outfmt 6
status=$?
echo -n " DONE "  | tee -a ${LOG}; date '+%H:%M:%S %Y-%m-%d' |tee -a ${LOG}

echo -n "# Extracting candidate PolA sequences ..............." | tee -a ${LOG}
cut -f2 "${OUT}.btab" | sort -u > "${OUT}_candidate_pola.lookup"
status=$?
extract_reads -f "${OUT}_candidate_pola.lookup" < "${FASTA}" > "${OUT}_candidate_pola.fasta"
status=$?
echo -n " DONE "  | tee -a ${LOG}; date '+%H:%M:%S %Y-%m-%d' |tee -a ${LOG}

echo -n "# Running the 762 caller ............................" | tee -a ${LOG}
762_caller.pl -i ${OUT}_candidate_pola.fasta \
    -r 16 \
    -o ${OUT}.762 \
    --fast
status=$?
echo -n " DONE "  | tee -a ${LOG}; date '+%H:%M:%S %Y-%m-%d' |tee -a ${LOG}



exit $status
