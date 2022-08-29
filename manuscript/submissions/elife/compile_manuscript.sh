#!/bin/bash

markdown_files=(*_AMBA_legend.Rmd)
markdown_files+=(*_figure*.Rmd)
markdown_files+=(*_gene_enrichment.Rmd)

echo "Generating figures..."

for file in ${markdown_files[@]};
do
	Rscript extract_code.R --infile $file
	
	r_file="${file/".Rmd"/".R"}"
	
	echo "Executing ${r_file} ..."
	Rscript $r_file 
	
	rm $r_file
	
	if [ -f Rplots.pdf ]; 
	then
		rm Rplots.pdf
	fi
done

echo "Rendering manuscript..."

Rscript render_pdf.R --infile TranscriptomicSimilarity_elife.Rmd

