if (!params.sequence){
	println "ERROR: Specify input nucleotide multifasta file"
	exit 1
}
if (!params.gencode){
	println "ERROR: Specify gencode number (e.g. 1,2,3 etc.)"
	exit 1
}
if (!params.outgroup){
	println "ERROR: Specify outgroup name of Id in input fasta file"
	exit 1
}
if (!params.aligned){
	println "ERROR: Specify if sequences are aligned or not"
	exit 1
}

if (!params.all){params.all = "false"}
if (!params.syn4f){params.syn4f = "false"}
if (!params.nonsyn){params.nonsyn = "false"}
if (!params.use_macse){params.use_macse = "false"} 
if (!params.verbose){params.verbose = "false"} 
if (!params.internal){params.internal = "false"} 
if (!params.terminal){params.terminal = "false"} 
if (!params.branch_spectra){params.branch_spectra = "false"}
if (!params.save_exp_mutations){params.save_exp_mutations = "false"}
if (!params.exclude_cons_sites){params.exclude_cons_sites = "false"}
if (!params.uncertainty_coef){params.uncertainty_coef = "false"}
if (!params.njobs){params.njobs = "1"}
if (!params.required_nseqs){params.required_nseqs = 4}
THREADS = params.njobs

// TODO add default values for params below

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

query_nucleotide_fasta = file(params.sequence, type: 'any') 
Channel.value(params.gencode).into{g_396_gencode_g_410;g_396_gencode_g_411;g_396_gencode_g_422;g_396_gencode_g_423;g_396_gencode_g_433}
Channel.value(params.aligned).set{aligned_param}
Channel.value(params.outgroup).set{outgroup_param}
Channel.value(params.use_macse).set{use_macse_param}

if (!params.treefile || params.treefile == ''){
	Channel.value("NO_FILE").set{precalculated_tree}
} else {
	precalculated_tree = file(params.treefile, type: 'any')
} 

// text_files = Channel.fromPath( '/path/*.txt' ).ifEmpty( file('./default.txt') ) for optioanl input

// TODO write log file fith full params list


process nucleotide_fasta_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /species_mapping.txt$/) "$filename"
}

input:
 file query from query_nucleotide_fasta
 val outgrp from outgroup_param

output:
 file "sequences.fasta"  into g_428_multipleFasta_g_433
 file "species_mapping.txt" into g_428_outputFileTxt

"""
nseqs=`grep -c ">" $query`
if [ \$nseqs -lt $params.required_nseqs ]; then
	echo "ERROR: Pipeline requires at least ${params.required_nseqs} sequences, but received \${nseqs}" >&2
	exit 1
fi

noutgrps=`grep -E -c "$outgrp" $query`
if [ \$noutgrps -eq 0 ]; then
	echo "Cannot find outgroup header in the query fasta." >&2
	exit 1
elif [ \$noutgrps -gt 1 ]; then
	echo "Cannot automatically find single outgroup in the query fasta. Found \$noutgrps records with '$outgrp' substring" >&2
	exit 1
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 5 char_numbers.log | grep -Ec "[ACGTacgt]"` -ge 3 ] && [ `grep -Ec "[EFILPQU]" char_numbers.log` -eq 0 ]; then
	echo "All right" >&2
else
	echo "Query fasta must contain nucleotides" >&2
	exit 1
fi

multifasta_coding.py -a $query -g "$outgrp" -o sequences.fasta -m species_mapping.txt
if [ ! -f species_mapping.txt ]; then
	echo 'required for compatibility reasons' > species_mapping.txt
fi

# TODO drop dublicates
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
 val aligned from aligned_param
 val use_macse from use_macse_param

output:
 file "msa_nuc.fasta"  into g_433_multipleFasta_g_420, g_433_multipleFasta_g_421, g_433_multipleFasta_g_422, g_433_multipleFasta_g_424, g_433_multipleFasta_g_425
 file "seq_dd_AA.fa" into seq_dd_AA_for_QC
 file "*.csv" optional true

"""
if [ $aligned = true ]; then
	# TODO add verification of aligment and delete this useless argument!!!!
	echo "No need to align!" >&2
	cp $seqs msa_nuc.fasta
	echo ">>>>>>>>>>>>>>>" > seq_dd_AA.fa

elif [ $aligned = false ]; then
	if [ $use_macse = true ]; then
		echo "Use macse as aligner" >&2
		java -jar /opt/macse_v2.07.jar -prog alignSequences -gc_def $gencode \
			-out_AA aln_aa.fasta -out_NT aln.fasta -seq $seqs
		echo "Do quality control" >&2
		/opt/scripts_latest/macse2.pl aln.fasta msa_nuc.fasta
		cp aln_aa.fasta seq_dd_AA.fa
	else
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
	fi
else
	echo "ERROR: Inappropriate value for 'aligned' parameter; must be 'true' or 'false'" >&2
	exit 1
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
├── species_mapping.txt					# Encoded headers of sequences (v2 for different version of pipeline)
├── logs/
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
 file prectree from precalculated_tree
 file labels_mapping from g_428_outputFileTxt

output:
 set val("iqtree"), file("iqtree.nwk")  into g_409_tree_g_315
 file "*.log" optional true  into g_409_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
# if input contains precalculated_tree
if [ $prectree != "input.2" ]; then
	echo "TODO need to correctly replace headers if needed"
	if [ ! -f $labels_mapping ]; then 
		echo "Something went wrong. Expected that tree must be relabeled \
		according to relabeled alignment, but file with labels map doesn't exist"
		exit 1
	fi

	awk '{print \$2 "\t" \$1}' $labels_mapping > species2label.txt
	nw_rename $prectree species2label.txt > iqtree.nwk

else
	iqtree2 -s $mulal -m $params.iqtree_model -nt $THREADS --prefix phylo
	mv phylo.treefile iqtree.nwk
	mv phylo.iqtree iqtree_report.log
	mv phylo.log iqtree.log
fi
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
