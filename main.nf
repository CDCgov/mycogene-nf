#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input    = ''
params.cyp51    = false
params.outdir   = "results"
params.query_aa = ''
params.query_fa = ''
params.help     = false

///// HELP MESSAGE /////
if (params.help) {
    help = """
|Usage:
|mycogene.nf --input <input_samplesheet> --outdir <output_dir>
|
|Required arguments:
| --input      Samplesheet with the Format: sampleID:path/to/fastq1:path/to/fastq2
| --outdir     Directory where process outputs are saved
| --query_aa   Amino acid sequence for the gene of interest
| --platform   Illumina/ONT
|
|Optional arguments:
| --cyp51      Run the analysis to identify SNPs and TR regions in the Cyp51 gene of all samples
| --query_fa   Nucleotide sequence for the gene of interest (required when --cyp51 is set)
| --help       Print this message and exit""".stripMargin()
    println(help)
    exit(0)
}

// Quality Filtering

process FASTP {
    tag "${sampleID}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy',
        saveAs: {filename -> filename.endsWith('.fastq.gz') ? filename : null}

    input:
    tuple val(sampleID), path(read1), path(read2)

    output:
    tuple val(sampleID), path("${sampleID}.R1.fastp.fastq.gz"), path("${sampleID}.R2.fastp.fastq.gz"), emit: trimmed
    tuple val(sampleID), path("${sampleID}_qc.json"),                                                   emit: json

    script:
    """
    fastp \
        -i ${read1} -I ${read2} \
        -o ${sampleID}.R1.fastp.fastq.gz \
        -O ${sampleID}.R2.fastp.fastq.gz \
        -e 20 \
        -j ${sampleID}_qc.json \
        -h /dev/null \
        --thread 4
    """
}

process FASTPLONG {
    tag "${sampleID}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy',
        saveAs: {filename -> filename.endsWith('.fastq.gz') ? "${sampleID}/${filename}" : null}

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.filtered.fastq.gz"), emit: trimmed
    tuple val(sampleID), path("${sampleID}_qc.json"),           emit: json

    script:
    """
    cat ${reads.join(' ')} > ${sampleID}.combined.fastq.gz
    fastplong \
        -i ${sampleID}.combined.fastq.gz \
        -o ${sampleID}.filtered.fastq.gz \
        -e 20 \
        -j ${sampleID}_qc.json \
        -h /dev/null
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
    tuple val(sampleID), path("${sampleID}/scaffolds.fasta"), emit: assembly

    script:
    """
    spades.py -1 ${read1} -2 ${read2} -k 127 --only-assembler -o ${sampleID}
    """
}

process ASSEMBLY_ONT {
    tag "${sampleID}"
    publishDir "${params.outdir}/Assemblies", mode: 'copy'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}/assembly.fasta"), emit: assembly

    script:
    """
    flye --nano-hq ${reads} -o ${sampleID}
    """
}

process blast_run {
    tag "${sampleID}"
    publishDir "${params.outdir}/blast_gene_hits", mode: 'copy', saveAs: { fname -> "${sampleID}/${fname}" }

    input:
    tuple val(sampleID), path(assembly)
    path(query_aa)

    output:
    tuple val(sampleID), path("prot_seq_${sampleID}.fasta"),    emit: prot_seq
    tuple val(sampleID), path("tblastn_raw_${sampleID}.tsv"),   emit: tblastn_raw

    script:
    """
    tblastn \
        -query ${query_aa} \
        -subject ${assembly} \
        -max_target_seqs 1 \
        -outfmt "6 delim=, qstart qend pident length sseq" \
    > tblastn_combined_${sampleID}.csv

    # Save qstart, qend, pident, length for coverage calculation
    sed 's/,/\t/g' tblastn_combined_${sampleID}.csv | cut -f1-4 > tblastn_raw_${sampleID}.tsv

    # Build the merged protein sequence (HSP overlap-aware)
    sed 's/,/\t/g' tblastn_combined_${sampleID}.csv \
    | sort -k1,1n \
    | awk '
        BEGIN { prev_end = 0 }
        {
            qstart = \$1; qend = \$2; seq = \$5
            if (qstart > prev_end) {
                printf "%s", seq
                prev_end = qend
            } else if (qend > prev_end) {
                printf "%s", substr(seq, prev_end - qstart + 2)
                prev_end = qend
            }
        }
        END { print "" }
    ' \
    | awk -v s="${sampleID}" 'BEGIN{print ">"s}{print}' \
    > prot_seq_${sampleID}.fasta
    """
}

// Gene coverage: identity-weighted coverage of the reference protein

process gene_coverage {
    tag "${sampleID}"

    input:
    tuple val(sampleID), path(tblastn_raw)
    val(ref_length)

    output:
    tuple val(sampleID), path("coverage_${sampleID}.csv"), emit: coverage

    script:
    """
    gene_coverage.py \
        --sample      ${sampleID} \
        --tblastn-out ${tblastn_raw} \
        --ref-length  ${ref_length} \
        --output      coverage_${sampleID}.csv
    """
}

// Per-sample QC report

process SAMPLE_QC {
    tag "${sampleID}"

    input:
    tuple val(sampleID), val(platform), path(json), path(assembly), path(coverage_csv)

    output:
    tuple val(sampleID), path("${sampleID}_qc.csv"), emit: qc_csv

    script:
    """
    sample_qc.py \
        --sample        ${sampleID} \
        --platform      ${platform} \
        --json          ${json} \
        --assembly      ${assembly} \
        --gene-coverage ${coverage_csv} \
        --output        ${sampleID}_qc.csv
    """
}

// Merge all per-sample QC CSVs into one run-level report

process MERGE_QC {
    publishDir "${params.outdir}/qc_report", mode: 'copy'

    input:
    path(qc_csvs)

    output:
    path("qc_report.csv")

    script:
    """
    merge_qc.py --input ${qc_csvs} --output qc_report.csv
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

process parse_mutations {
    publishDir "${params.outdir}/mutation_report", mode: 'copy'

    input:
    path(alignment)

    output:
    path("mutations.csv")

    script:
    """
    parse_mutations.py --alignment ${alignment} --output mutations.csv
    """
}

/// Cyp51 Analysis ONLY

process extract_best_hit {
    tag "${sampleID}"
    publishDir "${params.outdir}/intermediate_outputs", mode: 'copy'

    input:
    tuple val(sampleID), path(assembly)
    path(query_fa)

    output:
    tuple val(sampleID), path("best_hit_${sampleID}.tsv"), emit: best_hit

    script:
    """
    blastn -query ${query_fa} -subject ${assembly} -outfmt "6 qseqid sseqid sstart send pident length evalue bitscore" | sed 's,^,'"${sampleID}"'\t,' | head -n 1 > best_hit_${sampleID}.txt
    awk 'BEGIN {OFS="\t"; print "sample","qseqid","sseqid","sstart","send","pident","length","evalue","bitscore"}' > best_hit_${sampleID}.tsv
    awk 'BEGIN{OFS="\t"} {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9}' best_hit_${sampleID}.txt >> best_hit_${sampleID}.tsv
    """
}

process extract_cyp51_coding_noncoding_sequence {
    tag "${sampleID}"
    publishDir "${params.outdir}/intermediate_outputs", mode: 'copy'

    input:
    tuple val(sampleID), path(assembly), path(blast_out)

    output:
    tuple val(sampleID), path("seq_${sampleID}.fasta"),        emit: coding_seq
    tuple val(sampleID), path("wnoncoding_${sampleID}.fasta"), emit: non_coding_seq

    script:
    """
    read scaff start end <<< "\$(awk 'NR==2{print \$3, \$4, \$5}' ${blast_out})"

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
            0)  tr="no_TR"         ;;
            34) tr="TR34"          ;;
            46) tr="TR46"          ;;
            53) tr="TR53"          ;;
            *)  tr="unknown_indel" ;;
        esac
        echo -e "\$sampleID\t\$tr" >> TR_report.tsv
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
        platform : "Please specify the sequencing platform: 'illumina' or 'ont'"
    ]
    def missing = required.findAll { k, msg -> !params[k] }
    if ( missing ) {
        error """
ERROR: Missing required parameter(s):
  ${missing.collect { k, msg -> "--${k}: ${msg}" }.join('\n  ')}
"""
    }

    if ( params.cyp51 && !params.query_fa ) {
        error "Missing required parameter for Cyp51 Analysis: --query_fa"
    }

    def query_aa_file = file(params.query_aa)
    def query_fa_file = params.query_fa ? file(params.query_fa) : null

    // Reference protein length, used as the coverage denominator
    def ref_aa_length = 0
    query_aa_file.eachLine { line ->
        if (!line.startsWith(">")) {
            ref_aa_length += line.trim().length()
        }
    }

    if (params.platform == 'illumina') {
        ch_samples = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
        qc_reads  = FASTP(ch_samples)
        assemblies = ASSEMBLY_ILLUMINA(qc_reads.trimmed)
    }

    if (params.platform == 'ont') {
        ch_samples = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, file("${row.folder}/*.fastq.gz")) }
        qc_reads  = FASTPLONG(ch_samples)
        assemblies = ASSEMBLY_ONT(qc_reads.trimmed)
    }

    blast_ch    = blast_run(assemblies.assembly, query_aa_file)
    coverage_ch = gene_coverage(blast_ch.tblastn_raw, ref_aa_length)

    // QC report: join QC json + assembly + gene coverage per sample
    qc_input = qc_reads.json
        .join(assemblies.assembly)
        .join(coverage_ch.coverage)
        .map { sid, json, asm, cov -> tuple(sid, params.platform, json, asm, cov) }

    sample_qc_ch = SAMPLE_QC(qc_input)
    MERGE_QC(sample_qc_ch.qc_csv.map { sid, csv -> csv }.collect())

    // Filter: only QC-passing samples into alignment and downstream
    passing_prot_seq = sample_qc_ch.qc_csv
        .filter { sid, csv ->
            def lines  = csv.text.readLines()
            def header = lines[0].split(',')
            def vals   = lines[1].split(',')
            def row    = [header, vals].transpose().collectEntries()
            row.qc_status == 'PASS'
        }
        .map { sid, csv -> sid }
        .join(blast_ch.prot_seq)
        .map { sid, fasta -> fasta }

    // Filter: only QC-passing samples into cyp51
    passing_assemblies = sample_qc_ch.qc_csv
        .filter { sid, csv ->
            def lines  = csv.text.readLines()
            def header = lines[0].split(',')
            def vals   = lines[1].split(',')
            def row    = [header, vals].transpose().collectEntries()
            row.qc_status == 'PASS'
        }
        .map { sid, csv -> sid }
        .join(assemblies.assembly)
        .map { sid, asm -> tuple(sid, asm) }

    // Protein alignment, visualization and mutation report (QC-passing only)
    prot_in  = passing_prot_seq.collectFile(name: 'all_aligned_protein.fasta')
    prot_aln = combine_and_align(query_aa_file, prot_in)
    visualize_snps(prot_aln)
    parse_mutations(prot_aln)

    // Cyp51 analysis on QC-passing samples only
    if ( params.cyp51 ) {
        log.info "Running Cyp51/TR analysis"

        best_hit     = extract_best_hit(passing_assemblies, query_fa_file)
        extract_ch   = passing_assemblies.join(best_hit.best_hit)
        cyp51_seq_ch = extract_cyp51_coding_noncoding_sequence(extract_ch)

        tr_ch   = cyp51_seq_ch.non_coding_seq
        dist_TR = identify_distance_TR(tr_ch).collect()

        combined_wnoncoding = tr_ch
            .map { sid, seq -> seq }
            .collectFile(name: 'wnoncoding_gene_multifasta.fasta')

        alignment = align_cyp51(query_fa_file, combined_wnoncoding)
        plot_TR(alignment)
        report_TR(dist_TR)
    }
}
