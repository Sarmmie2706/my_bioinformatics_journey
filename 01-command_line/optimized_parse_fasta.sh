#!/bin/bash

organism=""
while IFS= read -r line
do
	if [[ ${line} == ">"* ]]
	then
		organism=$(echo ${line} | grep -oP "(?<=OS=).*?(?= OX)" | sed "s/ /_/g")
		if [[ -e "${organism}.fasta" ]]
		then
			echo ${line} >> "${organism}.fasta"
		else
			echo ${line} > "${organism}.fasta"
		fi
	else
		echo ${line} >> "${organism}.fasta"
	fi
done < uniprot_accession.fasta
