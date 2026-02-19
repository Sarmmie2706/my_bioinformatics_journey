#!/bin/bash

set -eu
mkdir organisms
organism=""
while IFS= read -r line
do
	if [[ ${line} == ">"* ]]
	then
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
	else
		echo "${line}" >> "organisms/${organism}.fasta"
	fi
done < uniprot_accession.fasta
