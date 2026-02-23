#!/bin/bash

# Rename file with your file name
file=sequence.gb

# If the file exits, extract unique gene names from the file
if [[ -e ${file} ]]
then
	gene_list=$(grep -oP '(?<=/gene=").*?(?=")' ${file} | sort -u)

	# Count and display the number of unique genes in the files
	count=$(echo "${gene_list}" | wc -l)
	echo "========================================"

	# Display based on the number of unique genes found. This is just for the sake of grammar though, lol
	if (( ${count} == 1 ))
	then
		echo "Found ${count} unique gene in ${file}:"
	else
		echo "Found ${count} unique genes in ${file}:"
	fi

	# Display the contents of the gene list
	echo "${gene_list}"
	echo "========================================"

	# Save the contents to a file named based on the genbank file
	#THis is to prevent multiple extensions
	echo "${gene_list}" > "${file::-3}_genelist.txt"
	echo ""
	echo "Gene list saved to ${file::-3}_genelist.txt"
else 
	echo "File doesn't exist"
fi

