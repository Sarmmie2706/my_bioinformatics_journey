#!/bin/bash

exec 3> optimize_parse_debug.log
export BASH_XTRACEFD=3
set -eux

organism=""
mkdir organisms

# Read every line in the input file
while IFS= read -r line
do
	# If the line is a header, extract the organism name to create a file,
	# else append line into existing file
	if [[ ${line} == ">"* ]]
	then
		organism=$(echo ${line} | grep -oP "(?<=OS=).*?(?= OX)" | sed "s/ /_/g")
		echo ${line} >> "organisms/${organism}.fasta"
	else
		echo ${line} >> "organisms/${organism}.fasta"
	fi
done < uniprot_accession.fasta

# Change into the created directory
cd organisms/

# Read every file in the new directory 
for file in *.fasta
do
	while IFS= read -r file_line
	do
		# If the file is a header, extract the organism name 
		if [[ ${file_line} == ">"* ]]
		then
			org_name=$(echo ${file_line} | grep -oP "(?<=_).*?(?= )")
		fi
	done < ${file}
	# Count the number of headers in each file
	count=$(less ${file} | grep -c ">")
        echo "There are ${count} ${org_name} sequences in the ${file} file"
done 

cd ..
# In order not to run into a "File already exists" error in the mkdir command while rerunning the script
# If you want to keep the parsed files, you can uncomment the next line 
# rm -r organisms/
