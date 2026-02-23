#!/bin/bash

exec 3> debug_multi_genbank_parser.log
export BASH_XTRACEFD=3
set -eux

for file in *.gb
do
        gene_list=$(grep -oP '(?<=/gene=").*?(?=")' ${file} | sort -u)
        # Count and display the number of unique genes in the files
        count=$(echo "${gene_list}" | wc -l)
        echo "========================================"
        if (( ${count} == 1 ))
        then
                echo "Found ${count} unique gene in ${file}:"
        else
                echo "Found ${count} unique genes in ${file}:"
        fi

        # Display the contents of the gene list
        echo "${gene_list}"

        # Create a folder based on the file's name and save the contents to a file named based on the genbank file
	mkdir ${file::-3}
	cd ${file::-3}
        echo "${gene_list}" > "${file::-3}_genelist.txt"
        echo ""
        echo "Gene list saved to ${file::-3}_genelist.txt"
	echo "========================================"
	echo ""
	echo ""

	# Move outside to the parent directory so that other files can be looped through
	cd ..
	
	# In case you want to run the file again, so as not to run into "File already exists" error, uncomment the next line
	# rm -r ${file::-3}
done
