#!/bin/bash

organism=""
mkdir organisms
while IFS= read -r line
do
	if [[ ${line} == ">"* ]]
	then
		organism=$(echo ${line} | grep -oP "(?<=OS=).*?(?= OX)" | sed "s/ /_/g")
		if [[ ! -e "organisms/${organism}.fasta" ]]
		then
			org_name=$(echo ${line} | grep -oP "(?<=_).*?(?= )")
			echo ${line} > "organisms/${organism}.fasta"
		else
			echo ${line} >> "organisms/${organism}.fasta"
		fi
	else
		echo ${line} >> "${organism}.fasta"
	fi
done < uniprot_accession.fasta

cd organisms/
for file in *.fasta
do
	while IFS= read -r file_line
	do
		if [[ ${file_line} == ">"* ]]
		then
			org_name=$(echo ${file_line} | grep -oP "(?<=_).*?(?= )")
		fi
	done < ${file}
	count=$(less ${file} | grep -c ">")
        echo "There are ${count} ${org_name} sequences in ${file}"
done 

cd ..
rm -r organisms/
