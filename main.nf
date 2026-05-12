#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = ''
params.cyp51 = false
params.outdir = "results"
params.query_aa = ''
params.query_fa = ''
params.help = false

///// HELP MESSAGE /////

if (params.help) {
        help = """
              |Usage: 
              |mycogene.nf --input <input_samplesheet> --outdir <output_dir>
              |       
              |Required arguments:     
              |  --input     Samplesheet with the Format: sampleID:path/to/fastq1:path/to/fastq2  
              |  --outdir    Directory where process outputs are saved     
              |  --query_aa  Amino acid sequence for the gene of interest
              |  --query_fa  Nucleotide sequence for the gene of interest
              |  --platform  Illumina/ONT 
              |
              |Optional arguments:  
              |  --cyp51     Run the analysis to identify SNPs and TR regions in the Cyp51 gene of all samples
              |  --help      Print this message and exit""".stripMargin()

    println(help)
    exit(0)
}

// Quality Filtering

process FAQCS {

    tag "${sampleID}"
    
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'

    input:
    tuple val(sampleID), path(read1), path(read2)

    output:
    tuple val(sampleID), path("${sampleID}/${sampleID}.1.trimmed.fastq"), path("${sampleID}/${sampleID}.2.trimmed.fastq"), emit: trimmed

    script:
    """
    FaQCs -1 ${read1} -2 ${read2} -q 30 --prefix ${sampleID} -d ${sampleID}
    """
}

process FASTPLONG {
    tag "${sampleID}"
    
    publishDir "${params.outdir}/filtered_reads", mode: 'copy',
        saveAs: { filename -> "${sampleID}/${filename}" }

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.filtered.fastq.gz"), emit: trimmed


    script:
    """
    cat ${reads.join(' ')} > ${sampleID}.combined.fastq.gz
    fastplong -i ${sampleID}.combined.fastq.gz -o ${sampleID}.filtered.fastq.gz -q 30 
    """

}

// Assembly

process ASSEMBLY_ILLUMINA {

    tag "${sampleID}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/Assemblies", mode: 'copy'

    input:
    tuple val(sampleID), path(read1), path(read2)

    output:
    tuple val(sampleID), path("${sampleID}/scaffolds.fasta"), emit: fasta

    script:
    """
    spades.py -1 ${read1} -2 ${read2} -k 127  --only-assembler -o ${sampleID}
    """

}

process ASSEMBLY_ONT {

    tag "${sampleID}"


    publishDir "${params.outdir}/Assemblies", mode: 'copy'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}/assembly.fasta"), emit: fasta

    script:
    """
    flye --nano-hq ${reads} -o ${sampleID}
    """

}

process blast_run {

    tag "${sampleID}"

    publishDir "${params.outdir}/blast_gene_hits", mode: 'copy', saveAs: { fname -> "${sampleID}/${fname}" }

    input:
    tuple val(sampleID), path(fasta)
    path(query_aa)
    path(query_fa)

    output:
    tuple val(sampleID), path("best_hit_${sampleID}.tsv"), emit: best_hit
    tuple val(sampleID), path("prot_seq_${sampleID}.fasta"), emit: prot_seq

    script:
    """
    #blastn - run and clean output

    blastn -query ${query_fa} -subject ${fasta} -outfmt "6 qseqid sseqid sstart send pident length evalue bitscore" | sed 's,^,'"${sampleID}"'\t,' | head -n 1 > best_hit_${sampleID}.txt
    awk 'BEGIN {printf "%-8s\\t%-10s\\t%-35s\\t%-8s\\t%-8s\\t%-8s\\t%-8s\\t%-8s\\t%-8s\\n", "sample", "qseqid", "sseqid", "sstart", "send", "pident", "length", "evalue", "bitscore"}' > best_hit_${sampleID}.tsv
    awk '{printf "%-8s\\t%-10s\\t%-35s\\t%-8s\\t%-8s\\t%-8s\\t%-8s\\t%-8s\\t%-8s\\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9}' best_hit_${sampleID}.txt >> best_hit_${sampleID}.tsv
    tblastn -query ${query_aa} -subject ${fasta} -max_target_seqs 1 -outfmt "6 delim=, qstart sseq" | sed 's/,/\t/g' | sort -k 1 | cut -f2 | tr -d '\n' | sed 's/EYCFLNRQ//' | awk -v s="${sampleID}" 'BEGIN{print ">"s}{print}' > prot_seq_${sampleID}.fasta
    """

}

/// Protein alignment analysis

process combine_and_align {

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path(query_aa)
    path(prot_seq)

    output:
    path("aln_protein_output.fasta")

    script:
    """
    cat ${query_aa} ${prot_seq} > protein_output.fasta
    clustalo -i protein_output.fasta -o aln_protein_output.fasta
    """

}

process visualize_snps {

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path(alignment)

    output:
    path("protein_aln_snp.html")

    script:
    """
    mview -in fasta -html head -coloring mismatch -colormap red ${alignment} > protein_aln_snp.html
    """

}

/// Cyp51 Analysis ONLY

process extract_cyp51_coding_noncoding_sequence {

    // errorStrategy 'ignore'

    tag "${sampleID}"

    publishDir "${params.outdir}/intermediate_outputs", mode: 'copy'

    input:
    tuple val(sampleID), path(assembly), path(blast_out)

    output:
    tuple val(sampleID), path("seq_${sampleID}.fasta"), emit: coding_seq
    tuple val(sampleID), path("wnoncoding_${sampleID}.fasta"), emit: non_coding_seq

    script:
    """
    #extract contig name, start and end for each sample
    read scaff start end <<< "\$(awk 'NR==2{print \$3, \$4, \$5}' ${blast_out})"

    #extract 500 upstream nucleotides because the repeat regions are present 500bp upstream of the cyp51 gene
    #start position in backed by 500 bp (-500 bp) before gene starts; incase gene is reverse match, 500 bp are added (+500 bp)

    #extract gene sequences
    if [ \$start -lt \$end ]; then 
        samtools faidx ${assembly} \$scaff:\$start-\$end -o seq_${sampleID}.fasta
        samtools faidx ${assembly} \$scaff:\$((\$start-500))-\$end | sed "1s/.*/>${sampleID}/" > wnoncoding_${sampleID}.fasta    
    else
        samtools faidx -i ${assembly} \$scaff:\$end-\$start -o seq_${sampleID}.fasta
        samtools faidx -i ${assembly} \$scaff:\$end-\$((\$start+500)) | sed "1s/.*/>${sampleID}/" > wnoncoding_${sampleID}.fasta
    fi
    """

}

process identify_distance_TR {

    tag "${sampleID}"

    publishDir "${params.outdir}/intermediate_outputs", mode: 'copy'

    input:
    tuple val(sampleID), path(non_coding_seq)

    output:
    path("distance_${sampleID}.txt")

    script:
    """
    blastn -query <(echo -e ">left\nAGAGTTGTCTAGAATCACGCGGTCC\n>right\nGGATGTGTGCTGAGCCGAATGAAAGTTGCCTAATTACTAAGGTGTAGTTC") -subject ${non_coding_seq} \
    -outfmt '6 qseqid sseqid length sstart send' -task blastn-short | sort -k2,2 -k 1,1 -u | paste - - | awk '{print "${sampleID}", \$2, \$9-\$5-1}' > distance_${sampleID}.txt
    """

}

process report_TR {
    publishDir "${params.outdir}/Cyp51_analysis/TR_report", mode: 'copy'

    input:
    path(dist_files)

    output:
    path("TR_report.tsv")

    script:
    """
    echo -e "sampleID\tTR" > TR_report.tsv

    for f in ${dist_files}; do
        read sampleID contig dist <<< \$(awk 'NR==1{print \$1, \$2, \$3}' "\$f")

        case "\$dist" in 
            0) tr="no_TR" ;;
            34) tr="TR34" ;;
            46) tr="TR46" ;;
            53) tr="TR53" ;;
            *) tr="unknown_indel" ;;
        esac

        echo -e "\$sampleID\\t\$tr" >> TR_report.tsv
    done
    """
}

process align_cyp51 {

    input:
    path(cyp51_gene)
    path(non_coding_seq)

    output:
    path("aln_wnoncoding_gene_output.fasta")

    script:
    """
    cat ${cyp51_gene} ${non_coding_seq} > gene_output.fasta
    clustalo -i gene_output.fasta -o aln_wnoncoding_gene_output.fasta
    """

}

process plot_TR {

    publishDir "${params.outdir}/Cyp51_analysis/TR_plot", mode: 'copy'

    input:
    path(alignment)

    output:
    path("data_upstream-CPY51_aln.html")

    script:
    """
    mview -in fasta -html head -css on -coloring identity ${alignment} > data_upstream-CPY51_aln.html
    """

}

workflow {

    def required = [
        input    : "Please specify the input samplesheet",
        query_aa : "Please specify the path to the gene amino acid sequence FASTA file",
        query_fa : "Please specify the path to the gene nucleotide sequence FASTA file",
        platform : "Please specify the sequencing platform: 'illumina' or 'ont'"
    ]

    def missing = required.findAll { k, msg -> !params[k] }
    if( missing ){
        error """
        ERROR: Missing required parameter(s):

        ${missing.collect {k, msg -> "--${k}: ${msg}"}.join('\n        ')}
        """
    }

    def query_aa_file = file(params.query_aa)
    def query_fa_file = file(params.query_fa)

    if (params.platform == 'illumina') {
        ch_samples = Channel.fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }

        qc_reads = FAQCS(ch_samples)
        assemblies = ASSEMBLY_ILLUMINA(qc_reads.trimmed)
    }

    if (params.platform == 'ont') {
        ch_samples = Channel.fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, files("${row.folder}/*.fastq.gz")) }

        qc_reads = FASTPLONG(ch_samples)
        assemblies = ASSEMBLY_ONT(qc_reads.trimmed)
    }
    

    blast_ch = blast_run(assemblies.fasta,query_aa_file,query_fa_file)

    // protein clustal alignment and visualization
    prot_in = blast_ch.prot_seq
              .map { sid, fasta -> fasta }
              .collectFile(name: 'all_aligned_protein.fasta')

    prot_aln = combine_and_align(query_aa_file, prot_in)
    visualize_snps(prot_aln)

    // // everything here onwards is exclusively cyp51 analysis

    // if ( params.cyp51 ) {

    //     log.info "Running Cyp51/TR analysis"

    //     extract_ch = assemblies.join(blast_ch.best_hit)
    //     cyp51_seq_ch = extract_cyp51_coding_noncoding_sequence(extract_ch)

    //     // extract non coding sequence part, identify TR in sample and combine results into one report
    //     tr_ch = cyp51_seq_ch.non_coding_seq
    //     dist_TR = identify_distance_TR(tr_ch).collect()
        
    //     // combine all the fasta with 500 bp upstream for the visual alignment
    //     combined_wnoncoding = tr_ch.map { sid, seq -> seq}
    //                           .collectFile(name: 'wnoncoding_gene_multifasta.fasta')
    //     alignment = align_cyp51(query_fa_file,combined_wnoncoding)
    //     plot_TR(alignment)


    //     //create report for the TR identified
    //     report_TR(dist_TR)

    // }
}