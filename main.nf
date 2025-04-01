#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress; bamstats } from "./lib/ingress"
include { configure_igv } from "./lib/common"
include { process_references } from "./subworkflows/process_references"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
MINIMAP_ARGS_PRESETS = [
    "dna": "-ax lr:hq -y",
    "rna": "-ax splice -uf -y"
]

// Create an MMI index
process makeMMIndex {
    label "wfalignment"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "combined_refs.mmi"
    cpus params.threads
    memory {
        def ref_size = combined_refs.size()
        combined_refs.size() > 1e9 ? "16 GB" : "11 GB"
    }
    input:
        path combined_refs, stageAs: "combined_refs.fasta"
        val minimap_args
    output:
        path "combined_refs.mmi"
    script:
    """
    minimap2 -t $task.cpus $minimap_args -d combined_refs.mmi combined_refs.fasta
    """
}

// Check if an MMI file contains the same references as the FASTA reference file.
process checkReferences {
    label "wfalignment"
    cpus params.threads
    memory {
        def ref_size = combined_refs.size()
        combined_refs.size() > 1e9 ? "16 GB" : "11 GB"
    }
    input:
        path "combined_references.mmi"
        path "combined_refs.fasta.fai"
        path combined_refs, stageAs: "combined_references.fasta"
    output:
        val true
    script:
    """
    # Read MMI references and check if they are in the FASTA fai file.
    workflow-glue check_reference_index --mmi_file combined_references.mmi --fasta_fai combined_refs.fasta.fai
    """
}

process alignReads {
    label "wfalignment"
    cpus params.threads
    memory {
        combined_refs.size() > 1e9 ? "16 GB" : "11 GB"
    }
    input:
        tuple val(meta), path(input)
        path combined_refs
        val is_xam
        val minimap_args
    output:
        tuple val(meta), path("aligned.sorted.bam"), path("aligned.sorted.bam.bai")
    script:
        int sorting_threads = Math.min((task.cpus / 3) as int, 3)
        int mapping_threads = task.cpus - sorting_threads
        // the minimum for `params.threads` in the schema is `4` and we should have
        // positive values for both thread vars, but can't hurt to make extra sure
        sorting_threads = Math.max(1, sorting_threads)
        mapping_threads = Math.max(1, mapping_threads)
    """
    ${is_xam ? "samtools fastq -T '*' $input" : "cat $input"} \
    | minimap2 -t $mapping_threads $minimap_args --cap-kalloc 100m --cap-sw-mem 50m \
        $combined_refs - \
    | samtools sort --write-index -@ ${sorting_threads - 1} -o reads.bam - \
        -o aligned.sorted.bam##idx##aligned.sorted.bam.bai
    """
}

process addStepsColumn {
    // TODO: we don't need 200 windows for very short references; find heuristics for
    // determining window length / number for such cases
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input: path "lengths.tsv"
    output: path "lengths_with_steps.tsv"
    """
    #!/usr/bin/env python
    import pandas as pd
    all = pd.read_csv('lengths.tsv', sep='\\t',
        dtype={"name": str, "lengths": int}
        )
    # the number of depth windows and maximum depth window size are determined
    # in `params.wf`
    all["step"] = (all["lengths"] // $params.wf.num_depth_windows).clip(
        1, $params.wf.max_depth_window_size
    )
    all.to_csv('lengths_with_steps.tsv', index=False, header=False, sep='\\t')
    """
}

process readDepthPerRef {
    // TODO: check if parallelisation with `xargs` or `parallel` is more efficient
    label "wfalignment"
    cpus 3
    memory "7 GB"
    input:
        tuple val(meta), path(alignment), path(index)
        path ref_len
    output:
        tuple val(meta), path("depth.all_regions.bed.gz"), env(max_depth_and_locus)
    script:
    """
    while IFS=\$'\\t' read -r name lengths steps; do
        mosdepth -n --fast-mode --by "\$steps" --chrom "\$name" -t $task.cpus \
            depth."\$name".temp $alignment \
        || echo "No alignments for "\$name""
        [[ -f depth."\$name".temp.regions.bed.gz ]] && \
            cat depth."\$name".temp.regions.bed.gz >> depth.all_regions.bed.gz
    done < $ref_len

    # remove all the temp files
    find -name 'depth.*.temp*' -delete

    # find the window with the highest depth (to be used as initial region to display in
    # the IGV panel)
    max_depth_and_locus=\$(
        workflow-glue get_max_depth_locus depth.all_regions.bed.gz \
            $params.wf.max_depth_window_size
    )
    """
}

process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-alignment-report.html"
    cpus 1
    memory {11.GB * task.attempt}
    maxRetries 1
    errorStrategy = 'retry'
    input:
        path "per-sample-data/*"
        path "refnames/*"
        path counts
        path "versions.txt"
        path "params.json"
    output:
        path "wf-alignment-report.html"
    script:
    String counts_args = (counts.name == OPTIONAL_FILE.name) ? "" : "--counts $counts"
    """
    workflow-glue report \
        --data per-sample-data \
        --refnames_dir refnames \
        --versions versions.txt \
        --params params.json \
        --wf-version $workflow.manifest.version \
        $counts_args
    """
}

process getVersions {
    label "wfalignment"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    cpus 1
    memory "2 GB"
    output:
        path "versions.txt"
    script:
    """
    python --version | tr -s ' ' ',' | tr '[:upper:]' '[:lower:]' > versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | (head -n 1 && exit 0) | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    ezcharts --version | sed 's/ /,/' >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    bgzip --version | head -n1 | sed -E 's/(.*) /\\1,/' >> versions.txt
    """
}

process getParams {
    label "wfalignment"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "params.json"
    cpus 1
    memory "2 GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process collectFilesInDir {
    label "wfalignment"
    cpus 1
    memory "2 GB"
    input: tuple val(dirname), path("staging_dir/*")
    output: path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}


// workflow module
workflow pipeline {
    take:
        sample_data
        refs
        counts
        depth_coverage
    main:
        // get params & versions
        workflow_params = getParams()
        software_versions = getVersions()

        // minimap2 args
        String minimap_args
        minimap_args = params.minimap_args ?: \
            MINIMAP_ARGS_PRESETS[params.minimap_preset]

        // handle references
        // if params.references contains MMI index file
        // use this as reference
        combined_mmi_file = Channel.of(OPTIONAL_FILE)
        // Process references although input is an MMI index
        // as Jbrowse needs the processed FASTA file
        refs = process_references(params.references)
        if (params.reference_mmi_file) {
            log.info("Using the provided MMI index as reference.")
            log.info("Indexing parameters (-k, -w or -H) will be overridden by parameters used in the prebuilt index.")
            minimap_reference = Channel.fromPath(params.reference_mmi_file, checkIfExists: true).first()
            // make sure mmi index contains the same references as the fasta
            checkReferences(minimap_reference, refs.combined_index, refs.combined)
        } else {
            minimap_reference = makeMMIndex(refs.combined, minimap_args)
        }

        if (params.bam) {
            ch_branched = sample_data.branch { meta, bam, bai, stats ->
                to_align: meta["is_unaligned"]
                aligned: true
            }
            ch_to_align = ch_branched.to_align
            | map { meta, bam, bai, stats -> [meta, bam] }
            // we ran ingress with `stats: false`, so we can drop `stats` (which would
            // only be `null`) here
            bam = ch_branched.aligned
            | map { meta, bam, bai, stats -> [meta, bam, bai] }
        } else {
            // FASTQ input
            ch_to_align = sample_data
            | map { meta, fastq, stats -> [meta, fastq] }
            bam = Channel.empty()
        }

        // run minimap
        bam = bam
        | mix(
            alignReads(ch_to_align, minimap_reference, params.bam as boolean, minimap_args)
        )

        // get stats
        stats = bamstats(bam, ["per_read_stats": params.per_read_stats])
        | multiMap { meta, bam, bai, stats ->
            // histograms and flagstat should always be there
            ArrayList hists = file(stats.resolve("*.hist"))
            Path flagstat = file(stats.resolve("bamstats.flagstat.tsv"))
            // readstats is only there when they were requested
            Path readstats = file(stats.resolve("bamstats.readstats.tsv.gz"))
            Path runids = file(stats.resolve("run_ids"))
            // keys must be defined after defs for ... reasons
            hists: [meta, hists]
            flagstat: [meta, flagstat]
            readstats: [meta, readstats.exists() ? readstats : null]
            runids: [meta, runids]
        }
        ArrayList ingressed_run_ids = []
        stats.runids.map{ meta, runids_f -> runids_f }.splitText().subscribe(
            onNext: {
                ingressed_run_ids += it.strip()
            },
            onComplete: {
                params.wf["ingress.run_ids"] = ingressed_run_ids
            }
        )

        // create a channel with the stats files needed for the report
        files_for_report = stats.hists
        | join(stats.flagstat)

        // get the sample names and sort them (will be used for IGV config below)
        sample_names = bam
        | map { meta, bam, bai -> meta.alias }
        | toSortedList

        // determine read_depth per reference / bam file if requested
        igv_locus = Channel.of(null)
        if (depth_coverage) {
            // add step column to ref lengths
            ref_lengths_with_steps = addStepsColumn(refs.lengths_combined)

            // get the depths
            readDepthPerRef(bam, ref_lengths_with_steps)
            | map { meta, depths, max_depth_and_locus -> [meta, depths] }

            // add to the channel with the files for the report
            files_for_report = files_for_report
            | join(
                readDepthPerRef.out
                | map { meta, depths, max_depth_and_locus -> [meta, depths] }
            )
            if (params.igv) {
                // decide which locus to show in the initial view of the IGV panel as
                // follows:
                // * get the locus with the largest depth for each sample
                // * iterate over the samples and pick the first locus that has more depth
                //   than `params.wf.igv_locus_depth_threshold`
                igv_locus = readDepthPerRef.out
                | map { meta, depths, max_depth_and_locus ->
                    (depth, locus) = max_depth_and_locus.split()
                    [meta.alias, depth as float, locus]
                }
                | toList
                // collect the max depths and loci in a map with the sample names as keys
                | map {
                    it.collectEntries { alias, depth, locus -> [alias, [locus, depth]] }
                }
                // combine the map with the list of sample names (as these are in the right
                // order); we use `cross` with an empty closure here because `