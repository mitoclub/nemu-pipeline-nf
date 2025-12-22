
NSEQS_LIMIT=33000
max_target_seqs = params.max_target_seqs

process tblastn_and_seqs_extraction {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sampled_sequences.fasta$/) "$filename"
	else if (filename =~ /report.blast$/) "logs/$filename"
	else if (filename =~ /.blast_output_*_filtered.tsv$/) "logs/$filename"
	else if (filename =~ /headers_mapping.txt$/) "$filename"
	else if (filename =~ /encoded_headers.txt$/) "$filename"
	else if (filename =~ /.*.taxids$/) "logs/$filename"
}

input:
 file query from g_398_multipleFasta_g_406
 val species_name from g_1_species_name_g_415
 val genus_taxid from genus_taxid_value
 val gencode from g_220_gencode_g_406
 val DB from g_15_commondb_path_g_406

output:
 file "sampled_sequences.fasta" into g_415_multipleFasta_g_409
 file "report.blast" optional true
 file "blast_output_*_filtered.tsv" optional true
 file "headers_mapping.txt" optional true
 file "encoded_headers.txt" optional true
 file "*.taxids" optional true

script:
"""
outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe"

echo "INFO: Collecting relatives taxids" >&2
if [[ "$species_name" == Homo* ]]; then
	get_species_taxids.sh -t 9604 > genus.taxids
else
	# TODO parse prepared table instead of fetching remote database. If table exist of course, make it general
	get_species_taxids.sh -t $genus_taxid  > genus.taxids
fi
if [ `cat genus.taxids | wc -l` -eq 0 ]; then
	echo "Internal server error during fetching of genus taxids. Try again later" >&2
	exit 1
fi
sleep 2

echo "INFO: Collecting species taxid information" >&2
get_species_taxids.sh -n "$species_name" | head -n 5 > sp_tax.info
if [ `cat sp_tax.info | wc -l` -eq 0 ]; then
	echo "Internal server error during fetching of species taxa information. Try again later" >&2
	exit 1
fi

if [[ `grep -e "rank : species" -e "rank : isolate" sp_tax.info` ]]; then 
	raw_sp_taxid=`grep Taxid sp_tax.info`
	species_taxid="\${raw_sp_taxid#*Taxid : }"
else
	echo "Cannot find species taxon id. Try to use another species name/taxid."  >&2
	exit 1
fi
sleep 1

echo "INFO: Collecting under-species taxids" >&2
get_species_taxids.sh -t \$species_taxid > species.taxids
if [ `cat species.taxids | wc -l` -eq 0 ]; then
	echo "Internal server error during fetching of species taxids. Try again later" >&2
	exit 1
fi

grep -v -f species.taxids genus.taxids > relatives.taxids

echo "INFO: Blasting species sequences in the nt" >&2
tblastn -db $DB -db_gencode $gencode -max_target_seqs $max_target_seqs \
		-query $query -out blast_output_species.tsv -evalue 0.00001 \
		-num_threads $THREADS -taxidlist species.taxids \
		-outfmt "\$outfmt"

echo "INFO: Filtering out bad hits: ident <= 80, query coverage <= 0.5" >&2
awk '\$2 > 80 && \$3 / \$4 > 0.5' blast_output_species.tsv > blast_output_species_filtered.tsv

echo "INFO: Checking required number of hits" >&2
nhits=`cat blast_output_species_filtered.tsv | wc -l`
if [ \$nhits -lt $params.required_nseqs ]; then
	echo "ERROR: there are only \$nhits valuable hits in the Nucleotide collection for given query," >&2
	echo "but needed at least ${params.required_nseqs}." >&2
	exit 1
fi

echo "INFO: Preparing coords for nucleotide sequences extraction" >&2
# entry|range|strand, e.g. ID09 x-y minus
awk '{print \$1, \$6, \$7, (\$NF ~ /^-/) ? "minus" : "plus"}' blast_output_species_filtered.tsv > raw_coords.txt
awk '\$2 > \$3 {print \$1, \$3 "-" \$2, \$4}' raw_coords.txt > coords.txt
awk '\$3 > \$2 {print \$1, \$2 "-" \$3, \$4}' raw_coords.txt >> coords.txt

echo "INFO: Getting species nucleotide sequences" >&2
blastdbcmd -db $DB -entry_batch coords.txt -outfmt %f -out species_sequences.fasta

echo "INFO: Checking required number of extracted seqs" >&2
nseqs=`grep -c '>' species_sequences.fasta`
if [ \$nseqs -lt $params.required_nseqs ]; then
	echo "ERROR: cannot extract more than \$nseqs seqs from the database for given query, but needed at least ${params.required_nseqs}" >&2
	exit 1
fi

NO_OUTGRP_MODE=
echo "INFO: Checking number of taxids for outgroup" >&2
if [ `cat relatives.taxids | wc -l` -eq 0 ]; then
	echo "WARNING: there are no taxids that can be used as outgroup." >&2
	echo "Maybe this species is single in the family, so pipeline can build incorrect phylogeneti tree" >&2
	# echo "due to potential incorrect tree rooting. You can select sequences and outgroup manually and" >&2
	# echo "run pipeline on your nucleotide sequences" >&2
	# exit 1
	NO_OUTGRP_MODE=1
else
	echo -e "INFO: Blasting for outgroup search" >&2
	tblastn -db $DB -db_gencode $gencode -max_target_seqs 10 \
			-query $query -out blast_output_genus.tsv -evalue 0.00001 \
			-num_threads $THREADS -taxidlist relatives.taxids \
			-outfmt "\$outfmt"
	
	echo "INFO: Filtering out bad hits: ident <= 70, query coverage <= 0.5" >&2
	awk '\$2 > 70 && \$3 / \$4 > 0.5' blast_output_genus.tsv | sort -rk 9 > blast_output_genus_filtered.tsv

	echo "INFO: Checking required number of hits for outgroup" >&2
	
	if [ `cat blast_output_genus_filtered.tsv | wc -l` -eq 0 ]; then
		echo "WARNING: there are no hits in the database that could be used as outgroup." >&2
		# echo "Unfortunately pipeline cannot analyse this species using nt databse." >&2
		# exit 1
		NO_OUTGRP_MODE=1
	else
		echo "INFO: Preparing outgroup coords for nucleotide sequences extraction" >&2
		# entry|range|strand, e.g. ID09 x-y minus
		head -n 1 blast_output_genus_filtered.tsv | awk '{print \$1, \$6, \$7, (\$NF ~ /^-/) ? "minus" : "plus"}' > raw_coords_genus.txt
		awk '\$2 > \$3 {print \$1, \$3 "-" \$2, \$4}' raw_coords_genus.txt >  coords_genus.txt
		awk '\$3 > \$2 {print \$1, \$2 "-" \$3, \$4}' raw_coords_genus.txt >> coords_genus.txt

		echo "INFO: Getting nucleotide sequences of potential outgroups" >&2
		blastdbcmd -db $DB -entry_batch coords_genus.txt -outfmt %f -out outgroup_sequence.fasta
	fi

fi

echo "INFO: Sequences headers encoding" >&2
if [ \$NO_OUTGRP_MODE ]; then
	echo "INFO: Preparing sequences without outgroup" >&2
	multifasta_coding.py -a species_sequences.fasta -g "string that cannot be outgroup" -o sampled_sequences.fasta -m encoded_headers.txt
else
	echo "INFO: Preparing sequences with outgroup" >&2
	ohead=`head -n 1 outgroup_sequence.fasta`
	cat outgroup_sequence.fasta species_sequences.fasta > merged.fasta
	multifasta_coding.py -a merged.fasta -g "\${ohead:1}" -o sampled_sequences.fasta -m encoded_headers.txt
fi	
"""
}


process duplicates_filtration {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /seqs_unique.fasta$/) "$filename"
}

input:
 file seqs from g_415_multipleFasta_g_409

output:
 file "seqs_unique.fasta"  into g_409_multipleFasta_g418_428

"""
/opt/scripts_latest/codon_alig_unique.pl $seqs 1>seqs_unique.fasta
"""
}


process nucleotide_fasta_qc {

input:
 file query from g_409_multipleFasta_g418_428

output:
 file "sequences.fasta"  into g_428_multipleFasta_g_433
 file "char_numbers.log"

"""
nseqs=`grep -c ">" $query`
if [ \$nseqs -lt $params.required_nseqs ]; then
	echo "ERROR: Pipeline requires at least ${params.required_nseqs} sequences, but received \${nseqs}" >&2
	exit 1
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 5 char_numbers.log | grep -Ec "[ACGTacgt]"` -ge 3 ] && [ `grep -Ec "[EFILPQU]" char_numbers.log` -eq 0 ]; then
	echo "All right" >&2
else
	echo "Query fasta must contain nucleotides" >&2
	exit 1
fi

mv $query sequences.fasta
"""
}
