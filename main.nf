#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input        = ''
params.add_sra_file = ''
params.cyp51        = false
params.outdir       = "results"
params.query_aa     = ''
params.query_fa     = ''
params.multi_query = ''
params.help         = false

///// HELP MESSAGE /////
if (params.help) {
    help = """
|Usage:
|mycogene.nf --input <input_samplesheet> --query_aa <gene_aa.fasta> --platform <illumina|ont> --outdir <output_dir>
|mycogene.nf --add_sra_file <accessions.txt> --query_aa <gene_aa.fasta> --platform <illumina|ont> --outdir <output_dir>
|mycogene.nf --input <input_samplesheet> --add_sra_file <accessions.txt> --query_aa <gene_aa.fasta> --platform <illumina|ont> --outdir <output_dir>
|mycogene.nf --input <input_samplesheet> --multi_query <genes.csv> --platform <illumina|ont> --outdir <output_dir>
|
|Every command needs a reads source, a query, --platform, and --outdir. See below.
|
|Provide at least one of:
| --input          Samplesheet (CSV, header row required):
|                   Illumina: sample,fastq_1,fastq_2
|                   ONT:      sample,folder  (folder containing that sample's *.fastq.gz)
| --add_sra_file   Text file of SRA run accessions (one SRR/ERR/DRR per line); reads are
|                   downloaded with sra-tools and combined with any --input samples
|
|Also required:
| --outdir     Directory where process outputs are saved
| --platform   Illumina/ONT
|
|Provide exactly one of:
| --query_aa      Amino acid sequence for a single gene of interest
| --multi_query   CSV file for multi-gene mode, no header, one gene per line:
|                  gene_name,/path/to/gene_protein.fasta
|                  Every sample is BLASTed against every gene in one run. Outputs keep
|                  the usual top-level folders (blast_gene_hits/, plots/, mutation_report/)
|                  with gene as a subfolder inside each, e.g. blast_gene_hits/<gene>/...
|                  Mutation files are flat: mutation_report/<gene>_mutations.csv.
|                  QC report stays a single results/qc_report/qc_report.csv with a 'gene'
|                  column, since coverage is evaluated per gene.
|                  NOTE: --cyp51 cannot currently be combined with --multi_query.
|
|Optional arguments:
| --cyp51      Run the analysis to identify SNPs and TR regions in the Cyp51 gene of all samples
|               (single-gene mode only)
| --query_fa   Nucleotide sequence for the gene of interest (required when --cyp51 is set)
| --help       Print this message and exit""".stripMargin()
    println(help)
    exit(0)
}

// Quality Filtering

process FASTP {
    tag "${sampleID}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'

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

// SRA download (combined with --input, if both are given)

process DOWNLOAD_SRA {
    tag "${accession}"
    publishDir "${params.outdir}/raw_reads", mode: 'copy'
    errorStrategy 'ignore'

    input:
    val(accession)

    output:
    tuple val(accession), path("*.fastq.gz")

    script:
    def split_flag     = params.platform == 'illumina' ? '--split-files' : ''
    def expected_files = params.platform == 'illumina' ? 2 : 1
    """
    prefetch ${accession} -O .
    vdb-validate ${accession}/${accession}.sra

    fasterq-dump ${accession}/${accession}.sra ${split_flag} --threads ${task.cpus} -O .

    shopt -s nullglob
    fastq_files=(*.fastq)
    if [ \${#fastq_files[@]} -ne ${expected_files} ]; then
        echo "ERROR: expected ${expected_files} fastq file(s) for ${accession} (platform=${params.platform}), got \${#fastq_files[@]}" >&2
        exit 1
    fi

    gzip "\${fastq_files[@]}"
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
    tag "${geneName ? geneName + ':' : ''}${sampleID}"
    publishDir "${params.outdir}/blast_gene_hits${geneName ? '/' + geneName : ''}", mode: 'copy', saveAs: { fname -> "${sampleID}/${fname}" }

    input:
    tuple val(sampleID), path(assembly), val(geneName), path(query_aa)

    output:
    tuple val(geneName), val(sampleID), path("prot_seq_${sampleID}.fasta"),  emit: prot_seq
    tuple val(geneName), val(sampleID), path("tblastn_raw_${sampleID}.tsv"), emit: tblastn_raw

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
    tag "${geneName ? geneName + ':' : ''}${sampleID}"

    input:
    tuple val(geneName), val(sampleID), path(tblastn_raw), val(ref_length)

    output:
    tuple val(geneName), val(sampleID), path("coverage_${sampleID}.csv"), emit: coverage

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
    tag "${geneName ? geneName + ':' : ''}${sampleID}"

    input:
    tuple val(sampleID), val(platform), val(geneName), path(json), path(assembly), path(coverage_csv)

    output:
    tuple val(geneName), val(sampleID), path("${geneName ? geneName + '_' : ''}${sampleID}_qc.csv"), emit: qc_csv

    script:
    """
    sample_qc.py \
        --sample        ${sampleID} \
        --platform      ${platform} \
        --gene          "${geneName}" \
        --json          ${json} \
        --assembly      ${assembly} \
        --gene-coverage ${coverage_csv} \
        --output        ${geneName ? geneName + '_' : ''}${sampleID}_qc.csv
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
    tag "${geneName}"
    publishDir "${params.outdir}/plots${geneName ? '/' + geneName : ''}", mode: 'copy'

    input:
    tuple val(geneName), path(query_aa), path(prot_seqs)

    output:
    tuple val(geneName), path("aln_protein_output.fasta")

    script:
    """
    cat ${query_aa} ${prot_seqs.join(' ')} > protein_output.fasta
    clustalo -i protein_output.fasta -o aln_protein_output.fasta
    """
}

process visualize_snps {
    tag "${geneName}"
    publishDir "${params.outdir}/plots${geneName ? '/' + geneName : ''}", mode: 'copy'

    input:
    tuple val(geneName), path(alignment)

    output:
    path("protein_aln_snp.html")

    script:
    """
    mview -in fasta -html head -coloring mismatch -colormap red ${alignment} > protein_aln_snp.html
    """
}

process parse_mutations {
    tag "${geneName}"
    publishDir "${params.outdir}/mutation_report", mode: 'copy'

    input:
    tuple val(geneName), path(alignment)

    output:
    path("${geneName ? geneName + '_' : ''}mutations.csv")

    script:
    """
    parse_mutations.py --alignment ${alignment} --output ${geneName ? geneName + '_' : ''}mutations.csv
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
        platform : "Please specify the sequencing platform: 'illumina' or 'ont'"
    ]
    def missing = required.findAll { k, msg -> !params[k] }
    if ( missing ) {
        error """
ERROR: Missing required parameter(s):
  ${missing.collect { k, msg -> "--${k}: ${msg}" }.join('\n  ')}
"""
    }

    if ( !params.input && !params.add_sra_file ) {
        error "ERROR: Specify at least one of --input (samplesheet) or --add_sra_file (SRA accession list)"
    }

    if ( !params.query_aa && !params.multi_query ) {
        error "ERROR: Specify one of --query_aa (single gene) or --multi_query (multi-gene CSV)"
    }

    if ( params.query_aa && params.multi_query ) {
        error "ERROR: Specify only one of --query_aa or --multi_query, not both"
    }

    if ( params.cyp51 && params.multi_query ) {
        error "ERROR: --cyp51 cannot currently be combined with --multi_query; run the Cyp51 TR analysis as a separate single-gene invocation"
    }

    if ( params.cyp51 && !params.query_fa ) {
        error "Missing required parameter for Cyp51 Analysis: --query_fa"
    }

    def query_fa_file = params.query_fa ? file(params.query_fa) : null

    // Build the gene channel: (geneName, queryAaFile, refAaLength)
    // Single-gene mode uses geneName = '' so publishDir/filenames stay exactly
    // as they were before multi-gene support (no /'' subfolder, no ''_ prefix).
    def countAaLength = { f ->
        def len = 0
        f.eachLine { line -> if (!line.startsWith('>')) len += line.trim().length() }
        return len
    }

    if ( params.multi_query ) {
        ch_genes = Channel.fromPath(params.multi_query)
            .splitCsv()
            .map { row ->
                def geneName  = row[0].trim()
                def queryFile = file(row[1].trim())
                tuple(geneName, queryFile, countAaLength(queryFile))
            }
    } else {
        def query_aa_file = file(params.query_aa)
        ch_genes = Channel.of( tuple('', query_aa_file, countAaLength(query_aa_file)) )
    }

    if ( params.add_sra_file ) {
        ch_accessions = Channel.fromPath(params.add_sra_file)
            .splitText()
            .map { it.trim() }
            .filter { it }

        dl_reads = DOWNLOAD_SRA(ch_accessions)
    }

    if (params.platform == 'illumina') {
        ch_from_input = params.input
            ? Channel.fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
            : Channel.empty()
        ch_from_sra = params.add_sra_file
            ? dl_reads.map { id, reads -> tuple(id, reads[0], reads[1]) }
            : Channel.empty()
        ch_samples = ch_from_input.mix(ch_from_sra)
        qc_reads  = FASTP(ch_samples)
        assemblies = ASSEMBLY_ILLUMINA(qc_reads.trimmed)
    }

    if (params.platform == 'ont') {
        ch_from_input = params.input
            ? Channel.fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample, file("${row.folder}/*.fastq.gz")) }
            : Channel.empty()
        ch_from_sra = params.add_sra_file
            ? dl_reads.map { id, reads -> tuple(id, reads) }
            : Channel.empty()
        ch_samples = ch_from_input.mix(ch_from_sra)
        qc_reads  = FASTPLONG(ch_samples)
        assemblies = ASSEMBLY_ONT(qc_reads.trimmed)
    }

    blast_input = assemblies.assembly.combine(ch_genes)   // (sampleID, assembly, geneName, queryFile, refLen)
        .map { sid, asm, gene, query, len -> tuple(sid, asm, gene, query) }

    blast_ch = blast_run(blast_input)

    ch_gene_reflen = ch_genes.map { gene, query, len -> tuple(gene, len) }

    coverage_input = blast_ch.tblastn_raw
        .combine(ch_gene_reflen, by: 0)                    // (geneName, sampleID, tsv, refLen)

    coverage_ch = gene_coverage(coverage_input)

    // QC report: broadcast per-sample json/assembly across every gene's coverage row
    qc_input = coverage_ch.coverage
        .map { gene, sid, cov -> tuple(sid, gene, cov) }
        .combine(qc_reads.json, by: 0)
        .combine(assemblies.assembly, by: 0)
        .map { sid, gene, cov, json, asm -> tuple(sid, params.platform, gene, json, asm, cov) }

    sample_qc_ch = SAMPLE_QC(qc_input)
    MERGE_QC(sample_qc_ch.qc_csv.map { gene, sid, csv -> csv }.collect())

    // Filter: only QC-passing (sample, gene) pairs go into alignment and downstream
    passing_keys = sample_qc_ch.qc_csv
        .filter { gene, sid, csv ->
            def lines  = csv.text.readLines()
            def header = lines[0].split(',')
            def vals   = lines[1].split(',')
            def row    = [header, vals].transpose().collectEntries()
            row.qc_status == 'PASS'
        }
        .map { gene, sid, csv -> tuple(gene, sid) }

    passing_prot_seq = passing_keys
        .combine(blast_ch.prot_seq, by: [0, 1])
        .map { gene, sid, fasta -> tuple(gene, fasta) }

    // Filter: only QC-passing samples into cyp51 (single-gene mode only, geneName == '')
    passing_assemblies = passing_keys
        .map { gene, sid -> sid }
        .join(assemblies.assembly)
        .map { sid, asm -> tuple(sid, asm) }

    // Protein alignment, visualization and mutation report (per gene, QC-passing only)
    align_input = ch_genes
        .map { gene, query, len -> tuple(gene, query) }
        .join(passing_prot_seq.groupTuple())

    prot_aln = combine_and_align(align_input)
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
