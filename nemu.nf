if (!params.sequence){
	println "INPUT ERROR: Specify input nucleotide multifasta file"
	exit 1
}
if (!params.gencode){
	println "INPUT ERROR: Specify gencode number (e.g. 1,2,3 etc.)"
	exit 1
}
if (!params.species_name){
	println "INPUT ERROR: Specify species name"
	exit 1
}
if (!params.DB){
	println "INPUT ERROR: Specify nucleotide database"
	exit 1
}
if (!params.iqtree_model){
	println "INPUT ERROR: Specify substitution model 'iqtree_model'"
	exit 1
}
if (!params.iqtree_anc_model){
	println "INPUT ERROR: Specify substitution model 'iqtree_anc_model'"
	exit 1
}
// if DB is NT, check that genus taxid specified
if (params.DB.endsWith("nt")){
	if (!params.genus_taxid){
		println "INPUT ERROR: Specify genus taxid when run pipeline on NT database"
		exit 1
	}
}

if (!params.all){params.all = "false"}
if (!params.syn4f){params.syn4f = "false"}
if (!params.nonsyn){params.nonsyn = "false"}
if (!params.run_shrinking){params.run_shrinking = "true"} 
if (!params.quantile){params.quantile = "0.1"} 
if (!params.verbose){params.verbose = "false"} 
if (!params.internal){params.internal = "false"} 
if (!params.terminal){params.terminal = "false"} 
if (!params.branch_spectra){params.branch_spectra = "false"}
if (!params.exclude_cons_sites){params.exclude_cons_sites = "false"}
if (!params.use_probabilities){params.use_probabilities = "false"}
if (!params.save_exp_mutations){params.save_exp_mutations = "false"}
if (!params.uncertainty_coef){params.uncertainty_coef = "false"}
if (!params.njobs){params.njobs = "1"}
if (!params.required_nseqs){params.required_nseqs = 3}
THREADS = params.njobs

// TODO add specific params logs
if (params.verbose == 'true') {
	println ""
	println "PARAMETERS:"
	println "all: ${params.all}"
	println "syn: true"
	println "syn4f: ${params.syn4f}"
	println "non-syn: ${params.nonsyn}"
	println "Use macse aligner: ${params.use_macse}"
	println "Minimal number of mutations to save 192-component spectrum (mnum192): ${params.mnum192}"
	println "Minimal number of sequences (leaves in a tree ingroup) to run the pipeline: ${params.required_nseqs}"
	println "Use probabilities: ${params.use_probabilities}"
	if (params.use_probabilities == 'true'){
		println "Mutation probability cutoff: ${params.proba_cutoff}"
	}
	println "Run simulation: ${params.run_simulation}"
	if (params.run_simulation == 'true'){
		println "Replics in simulation: ${params.replics}"
		println "Tree scaling coefficient: ${params.scale_tree}"
	}
	println "Threads: ${THREADS}"
	println "Run treeshrink: ${params.run_shrinking}"
	if (params.run_shrinking == 'true'){
		println "Shrinking Quantile: ${params.quantile}"
	}
	println "Exclude conservative sites: ${params.exclude_cons_sites}"
	println "internal branches spectra: ${params.internal}"
	println "terminal branches spectra: ${params.terminal}"
	println "branch-scpecific spectra: ${params.branch_spectra}"
	println ""
}


params.all_arg   = params.all == "true" ? "--all" : ""
params.syn4f_arg = params.syn4f == "true" ? "--syn4f" : ""
params.nonsyn_arg = params.nonsyn == "true" ? "--nonsyn" : ""
params.proba_arg = params.use_probabilities == "true" ? "--proba" : ""

Channel.value(params.DB).into{g_15_commondb_path_g_406;g_15_commondb_path_g_444}
Channel.value(params.genus_taxid).set{genus_taxid_value}
query_protein_sequence = file(params.sequence, type: 'any') 
Channel.value(params.gencode).into{g_220_gencode_g_406;g_396_gencode_g_410;g_396_gencode_g_411;g_396_gencode_g_422;g_396_gencode_g_423;g_396_gencode_g_433}
Channel.value(params.species_name).set{g_1_species_name_g_415}


process query_qc {

input:
 file query from query_protein_sequence

output:
 file "query_single.fasta"  into g_398_multipleFasta_g_406
 file "input_seq_char_counts.log" optional true  into g_398_logFile

"""
echo "Init query QC" >&2
if [ `grep -c ">" $query` -ne 1 ]; then
	echo "ERROR: Query fasta must contain single amino acid sequence" >&2
	exit 1
else
	echo "INFO: Number of sequences: `grep -c '>' $query`" >&2
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > input_seq_char_counts.log
if [ `head -n 4 input_seq_char_counts.log | grep -Ec "[ACGT]"` -lt 4 ] || [ `grep -Ec "[EFILPQU]" input_seq_char_counts.log` -ne 0 ]; then
	echo "INFO: It's probably amino asid sequence" >&2
else
	echo "ERROR: Query fasta must contain single amino acid sequence" >&2
	exit 1
fi

mv $query query_single.fasta
"""
}


NSEQS_LIMIT=33000
max_target_seqs = 1000

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

if [[ $DB == *nt ]]; then
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
	if [[ `grep "rank : species" sp_tax.info` ]]; then 
		raw_sp_taxid=`grep Taxid sp_tax.info`
		species_taxid="\${raw_sp_taxid#*Taxid : }"
	fi
	sleep 1

	echo "INFO: Collecting under-species taxids" >&2
	get_species_taxids.sh -t \$species_taxid > species.taxids
	if [ `cat species.taxids | wc -l` -eq 0 ]; then
		echo "Internal server error during fetching of species taxids. Try again later" >&2
		exit 1
	fi

	grep -v -f species.taxids genus.taxids > relatives.taxids

	echo "INFO: Checking number of taxids for outgroup" >&2
	if [ `cat relatives.taxids | wc -l` -eq 0 ]; then
		echo "ERROR: there are no taxids that can be used as outgroup." >&2
		echo "Maybe this species is single in the family, so pipeline can build incorrect phylogeneti tree" >&2
		echo "due to potential incorrect tree rooting. You can select sequences and outgroup manually and" >&2
		echo "run pipeline on your nucleotide sequences" >&2
		exit 1
	fi

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

	echo "INFO: Getting nucleotide sequences" >&2
	blastdbcmd -db $DB -entry_batch coords.txt -outfmt %f -out nucleotide_sequences.fasta

	echo "INFO: Checking required number of extracted seqs" >&2
	nseqs=`grep -c '>' nucleotide_sequences.fasta`
	if [ \$nseqs -lt $params.required_nseqs ]; then
		echo "ERROR: cannot extract more than \$nseqs seqs from the database for given query, but needed at least ${params.required_nseqs}" >&2
		exit 1
	fi

	echo -e "INFO: Blasting for outgroup search" >&2
	tblastn -db $DB -db_gencode $gencode -max_target_seqs 10 \
			-query $query -out blast_output_genus.tsv -evalue 0.00001 \
			-num_threads $THREADS -taxidlist relatives.taxids \
			-outfmt "\$outfmt"
	
	echo "INFO: Filtering out bad hits: ident <= 70, query coverage <= 0.5" >&2
	awk '\$2 > 70 && \$3 / \$4 > 0.5' blast_output_genus.tsv | sort -rk 9 > blast_output_genus_filtered.tsv

	echo "INFO: Checking required number of hits for outgroup" >&2
	if [ `cat blast_output_genus_filtered.tsv | wc -l` -eq 0 ]; then
		echo "ERROR: there are no hits in the database that could be used as outgroup." >&2
		echo "Unfortunately pipeline cannot analyse this species using nt databse." >&2
		exit 1
	fi

	echo "INFO: Preparing outgroup coords for nucleotide sequences extraction" >&2
	# entry|range|strand, e.g. ID09 x-y minus
	head -n 1 blast_output_genus_filtered.tsv | awk '{print \$1, \$6, \$7, (\$NF ~ /^-/) ? "minus" : "plus"}' > raw_coords_genus.txt
	awk '\$2 > \$3 {print \$1, \$3 "-" \$2, \$4}' raw_coords_genus.txt >  coords_genus.txt
	awk '\$3 > \$2 {print \$1, \$2 "-" \$3, \$4}' raw_coords_genus.txt >> coords_genus.txt

	echo "INFO: Getting nucleotide sequences of potential outgroups" >&2
	blastdbcmd -db $DB -entry_batch coords_genus.txt -outfmt %f -out outgroup_sequence.fasta
	
	echo "INFO: Sequences headers encoding" >&2
	ohead=`head -n 1 outgroup_sequence.fasta`
	cat outgroup_sequence.fasta nucleotide_sequences.fasta > sample.fasta
	multifasta_coding.py -a sample.fasta -g "\${ohead:1}" -o sampled_sequences.fasta -m encoded_headers.txt

else
	report=report.blast
	nseqs=$max_target_seqs

	while true
	do   
		echo "INFO: Blasting in midori2 database; max_target_seqs=\$nseqs" >&2
		tblastn -db $DB -db_gencode $gencode -num_descriptions \$nseqs -num_alignments \$nseqs \
				-query $query -out \$report -num_threads $THREADS
		
		if [ `grep -c "No hits found" \$report` -eq 0 ]; then 
			echo "INFO: some hits found in the database for given query" >&2
		else
			echo "ERROR: there are no hits in the database for given query" >&2
			exit 1
		fi

		if [ `grep -c "$species_name" \$report` -ge \$((nseqs * 2 - 10)) ]; then
			nseqs=\$((nseqs * 4))
			if [ \$nseqs -gt $NSEQS_LIMIT ]; then
				echo "UNEXPECTED ERROR: database cannot contain more than $NSEQS_LIMIT sequences of one gene of some species" >&2
				exit 1
			fi
			echo "INFO: run blasting again due to abcence of different species; nseqs=\$nseqs" >&2
		else
			echo "SUCCESS: other species for outgroup selection is in hits" >&2
			break
		fi
	done

	mview -in blast -out fasta \$report 1>raw_sequences.fasta

	/opt/scripts_latest/header_sel_mod3.pl raw_sequences.fasta "$species_name" 1>useless_seqs.fasta 2>headers_mapping.txt

	/opt/scripts_latest/nuc_coding_mod.pl headers_mapping.txt $DB 1>sampled_sequences.fasta

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

if [ `grep -E -c ">OUTGRP" $query` -ne 1 ]; then
	echo "Unexpected error. Cannot find outgroup header in the alignment." >&2
	echo "Probably outgroup sequence for this gene cannot be found in the database." >&2
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

thr_gaps = 0.05

process MSA {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /msa_nuc.fasta$/) "$filename"
	else if (filename =~ /.*.csv$/) "logs/$filename"
}

input:
 file seqs from g_428_multipleFasta_g_433
 val gencode from g_396_gencode_g_433
 val DB from g_15_commondb_path_g_444

output:
 file "msa_nuc.fasta"  into g_433_multipleFasta_g_420, g_433_multipleFasta_g_421, g_433_multipleFasta_g_422, g_433_multipleFasta_g_424, g_433_multipleFasta_g_425
 file "seq_dd_AA.fa" into seq_dd_AA_for_QC
 file "*.csv" optional true

"""
if [[ $DB == *nt ]]; then
	mafft --thread $THREADS $seqs > seqM.fa
	sed '/^>/!s/[actg]/\\U&/g' seqM.fa > seqMU.fa
	goalign clean seqs -c 0.3 -i seqMU.fa -o seqMC.fa
	goalign clean sites -c $thr_gaps -i seqMC.fa -o seqMCC.fa
	seqkit rmdup -s < seqMCC.fa > seq_dd.fa
	java -jar /opt/macse_v2.07.jar -prog alignSequences -seq seq_dd.fa \
		-gc_def $gencode -optim 2 -max_refine_iter 0 -ambi_OFF
	
	java -jar /opt/macse_v2.07.jar -prog exportAlignment \
		-align seq_dd_NT.fa -gc_def $gencode -ambi_OFF \
		-codonForInternalStop "NNN" -codonForFinalStop "---" \
		-codonForInternalFS "NNN" -codonForExternalFS "---" \
		-out_stat_per_seq macse_stat_per_seq.csv -out_stat_per_site macse_stat_per_site.csv 

	sed 's/!/n/g' seq_dd_NT_NT.fa > seq_dd_NT_FS.fa
	goalign clean sites -c $thr_gaps -i seq_dd_NT_FS.fa -o seq_dd_NT_FS_clean.fa
	seqkit rmdup -s < seq_dd_NT_FS_clean.fa > msa_nuc_lower.fasta
	sed '/^>/!s/[actg]/\\U&/g' msa_nuc_lower.fasta > msa_nuc.fasta

else
	# NT2AA
	java -jar /opt/macse_v2.07.jar -prog translateNT2AA -seq $seqs \
		-gc_def $gencode -out_AA translated.faa
	#ALN AA
	mafft --thread $THREADS translated.faa > translated_aln.faa
	#AA_ALN --> NT_ALN
	java -jar /opt/macse_v2.07.jar -prog reportGapsAA2NT \
		-align_AA translated_aln.faa -seq $seqs -out_NT aln.fasta
	echo "Do quality control" >&2
	/opt/scripts_latest/macse2.pl aln.fasta msa_nuc.fasta

	cp translated_aln.faa seq_dd_AA.fa
fi
"""
}


process MSA_QC {

input:
 file msa from g_433_multipleFasta_g_422
 file aa from seq_dd_AA_for_QC

"""
nstops=`grep -Eo "\\*[A-Za-z]" $aa | wc -l`
nseqs_aa=`grep -c ">" $aa`
thr_for_nstops=\$((nseqs_aa * 2)) # num of seqs * 2 is the max number of stops
if [ \$nstops -gt \$thr_for_nstops ]; then
	echo "There are stops in \${nstops} sequences. It's possible that you set incorrect gencode" >&2
	exit 1
fi

nseqs=`grep -c ">" $msa`
if [ \$nseqs -lt $params.required_nseqs ]; then
	echo "ERROR: Too low number of sequences! Pipeline requires at least ${params.required_nseqs} sequences, but after deduplication left only \${nseqs}" >&2
	exit 1
fi
"""
}


process write_readme {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /readme.txt$/) "$filename"
}

input:
 file something from g_433_multipleFasta_g_425
output:
 file "readme.txt"

"""
cat > readme.txt <<- EOM
Output structure:

.
├── final_tree.nwk						# Final phylogenetic tree
├── seqs_unique.fasta					# Filtered orthologous sequences
├── msa_nuc.fasta						# Verified multiple sequence alignment
├── headers_mapping.txt					# Encoded headers of sequences
├── encoded_headers.txt					# Encoded headers of sequences (v2 for different versions of input)
├── logs/
│   ├── report.blast					# Tblastn output during orthologs search
│   ├── *.taxids						# Taxids used in taxa-specific blasing in nt; relatives.taxids contains 
│	│									# 	other species from the genus of query and used for outgroup selection
│   ├── iqtree.log						# IQ-TREE logs during phylogenetic tree inference
│   ├── iqtree_report.log				# IQ-TREE report during phylogenetic tree inference
│   ├── iqtree_treeshrink.log			# TreeShrink logs
│   ├── iqtree_pruned_nodes.log			# Nodes pruned from tree by TreeShrink
│   ├── iqtree_anc.log					# IQ-TREE logs during ancestral reconstrution
│   ├── iqtree_anc_report.log			# IQ-TREE report during ancestral reconstrution
│   ├── iqtree_mut_extraction.log		# Logs during mutation extraction process
│   └── branches.txt					# Tree branch lenghts
├── figures
│   ├── ms12syn.pdf						# Barplot with  12-component spectrum on synonymous mutations
│   └── ms192syn.pdf					# Barplot with 192-component spectrum on synonymous mutations
├── tables
│   ├── rates.tsv						# Site rates categories for an alignment
│   ├── expected_freqs.tsv				# Frequencies of substitutions for each tree node genome
│   ├── mean_expexted_mutations.tsv		# Averaged frequencies of substitutions for entire tree
│   ├── ms12syn.tsv						# table with 12-component spectrum on synonymous mutations
│   ├── ms192syn.tsv					# table with 192-component spectrum on synonymous mutations
│   └── observed_mutations.tsv			# Recontructed mutations
EOM
"""
}


process fasta2states_table {

input:
 file aln from g_433_multipleFasta_g_421

output:
 file "leaves_states.state"  into g_421_state_g_410

"""
alignment2iqtree_states.py $aln leaves_states.state
"""
}


process ML_tree_inference {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file mulal from g_433_multipleFasta_g_420

output:
 set val("iqtree"), file("iqtree.nwk")  into g_409_tree_g_315
 file "*.log" optional true  into g_409_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
iqtree2 -s $mulal -m $params.iqtree_model -nt $THREADS --prefix phylo
mv phylo.treefile iqtree.nwk
mv phylo.iqtree iqtree_report.log
mv phylo.log iqtree.log
"""
}


process shrink_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g_409_tree_g_315

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_315_tree_g_302, g_315_tree_g_132
 file "*.log" optional true into g_315_logFile

"""
if [ $params.run_shrinking = true ] && [ `nw_stats $tree | grep leaves | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $params.quantile -x OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	if [ $params.run_shrinking = true ]; then
		echo "Shrinking are useless on such a small number of sequences" > ${name}_treeshrink.log
	fi
fi
"""
}


process terminal_branch_lengths_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /branches.txt$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g_315_tree_g_132

output:
 set val(name), file("branches.txt")  into g_132_branches

"""
nw_distance -m p -s f -n $tree | sort -grk 2 1> branches.txt

if [ `grep OUTGRP branches.txt | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` = False ]; then
	cat "branches.txt"
	echo "Something went wrong: outgroup is not furthest leaf in the tree" >&2
	exit 1
fi
"""
}


process tree_rooting_iqtree {

input:
 set val(name), file(tree) from g_315_tree_g_302

output:
 set val("${name}_rooted"), file("*.nwk")  into g_302_tree_g_326

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk
"""
}


estimate_rates = params.exclude_cons_sites == "true" ? "--rate" : ""

process ASR {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /final_tree.nwk$/) "$filename"
	else if (filename =~ /rates.tsv$/) "tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file mulal from g_433_multipleFasta_g_424
 set val(namet), file(tree) from g_302_tree_g_326

output:
 set val("iqtree"), file("final_tree.nwk")  into g_326_tree_g_410, g_326_tree_g_422
 set val("iqtree"), file("iqtree_anc.state")  into g_326_state_g_410
 path "rates.tsv" into g_326_ratefile
 file "*.log"  into g_326_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
nw_labels -I $tree | sed 's/\$/\$/' > leaves.txt
select_records.py -f leaves.txt --fmt fasta $mulal msa_filtered.fasta

iqtree2 -te $tree -s msa_filtered.fasta -m $params.iqtree_anc_model -asr -nt $THREADS --prefix anc $estimate_rates
if [ ! -f anc.rate ]; then
	touch anc.rate
fi
mv anc.rate rates.tsv
mv anc.iqtree iqtree_anc_report.log
mv anc.log iqtree_anc.log
nw_reroot anc.treefile OUTGRP | sed 's/;/ROOT;/' > final_tree.nwk

iqtree_states_add_part.py anc.state iqtree_anc.state
"""
}


save_exp_muts = params.save_exp_mutations == "true" ? "--save-exp-muts" : ""
use_uncertainty_coef = params.uncertainty_coef == "true" ? "--phylocoef" : "--no-phylocoef"

process mutations_reconstruction {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(namet), file(tree) from g_326_tree_g_410
 set val(label), file(states1) from g_326_state_g_410
 file states2 from g_421_state_g_410
 val gencode from g_396_gencode_g_410
 path rates from g_326_ratefile

output:
 file "observed_mutations.tsv"  into g_410_outputFileTSV
 file "expected_freqs.tsv"  into g_411_outputFileTSV
 file "expected_mutations.tsv" optional true
 file "mut_extraction.log"

"""
if [ $params.exclude_cons_sites = true ]; then 
	collect_mutations.py --tree $tree --states $states1 --states $states2 \
		--gencode $gencode --syn $params.syn4f_arg $params.nonsyn_arg \
		$params.proba_arg --no-mutspec --pcutoff $params.proba_cutoff \
		--mnum192 $params.mnum192 --outdir mout $save_exp_muts $use_uncertainty_coef \
		--rates $rates --cat-cutoff $params.cons_cat_cutoff
else
	collect_mutations.py --tree $tree --states $states1 --states $states2 \
		--gencode $gencode --syn $params.syn4f_arg $params.nonsyn_arg \
		$params.proba_arg --no-mutspec --pcutoff $params.proba_cutoff \
		--mnum192 $params.mnum192 --outdir mout $save_exp_muts $use_uncertainty_coef
fi
mv mout/* .
mv mutations.tsv observed_mutations.tsv
mv run.log mut_extraction.log
"""
}


process spectra_calculation {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.pdf$/) "figures/$filename"
}

input:
 file obs_muts from g_410_outputFileTSV
 file exp_freqs from g_411_outputFileTSV

output:
 file "*.tsv"
 file "*.pdf" optional true

"""
nmuts=`cat $obs_muts | wc -l`
if [ \$nmuts -lt 2 ]; then
	echo "ERROR: There are no reconstructed mutations after pipeline execution." >&2
	echo "Unfortunately this gene cannot be processed authomatically on available data." >&2
fi

calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
	--exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
	--proba_cutoff $params.proba_cutoff --plot -x pdf \
	--syn $params.syn4f_arg $params.all_arg $params.nonsyn_arg

if [ $params.internal = true ]; then
	calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
		--proba_cutoff $params.proba_cutoff --plot -x pdf --subset internal \
		--syn $params.syn4f_arg $params.all_arg $params.nonsyn_arg
	rm mean_expexted_mutations_internal.tsv
fi
if [ $params.terminal = true ]; then
	calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
		--proba_cutoff $params.proba_cutoff --plot -x pdf --subset terminal \
		--syn $params.syn4f_arg $params.all_arg $params.nonsyn_arg
	rm mean_expexted_mutations_terminal.tsv
fi
if [ $params.branch_spectra = true ]; then
	calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
		--proba_cutoff $params.proba_cutoff --branches \
		--syn $params.syn4f_arg $params.all_arg $params.nonsyn_arg
fi
"""
}


if (params.verbose == 'true') {
	workflow.onComplete {
	println "##Completed at: ${workflow.complete}; Duration: ${workflow.duration}; Success: ${workflow.success ? 'OK' : 'failed' }"
	}
}
