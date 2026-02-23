#!/bin/bash

# Adding file to capture executions
exec 3> parse_debug.log
export BASH_XTRACEFD=3
set -eux
mkdir organisms
organism=""

# Read every line in the input file
while IFS= read -r line
do
	# If the line is a header line, set the fasta file to input the line into
	# If it isn't, append to the current fasta file 
	if [[ ${line} == ">"* ]]
	then
		# Checking which organism it is to determine which fasta file is used
		case ${line} in
		*"OS=Pan"*)
		organism="chimp"
		;;
		*"OS=Homo"*)
                organism="human"
                ;;
		*"OS=Drosophila"*)
                organism="fruitfly"
                ;;
		*"OS=Mus"*)
                organism="mouse"
                ;;
		*"OS=Bos"*)
                organism="bovine"
                ;;
		*"OS=Sus"*)
                organism="pig"
                ;;
		*"OS=Rattus"*)
                organism="rat"
                ;;
		*"OS=Canis"*)
                organism="dog"
                ;;
		*"OS=Xenopus"*)
                organism="xenopus"
                ;;
		*)
		organism="animal"
		esac
		echo "${line}" >> "organisms/${organism}.fasta"
		echo "Done setting new header's fasta file"
	else
		echo "${line}" >> "organisms/${organism}.fasta"
	fi
done < uniprot_accession.fasta

# To prevent file a "File already exists error" if you want to rerun the script, uncomment the last line
# rm -r organisms
