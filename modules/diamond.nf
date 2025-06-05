process make_diamond_db {

    tag "${id}"

    input:
    tuple val( id ), path( fasta )

    output:
    tuple val( id ), path( "db.dmnd" )

    script:
    """
    diamond makedb --in "${fasta}" -d db.dmnd
    
    """

}


process diamond_blastp {

    tag "${id}"

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${id}-${it}" }
    )

    input:
    tuple val( id ), path( db ), path( queries )

    output:
    tuple val( id ), path( "hits.tsv" ), emit: data
    tuple val( id ), path( "hit_count.tsv" ), emit: stats

    script:
    """
    header="target_uniprot_id,species_target_uniprot_id,target_accession,species_target_accession,query_length,target_length,alignment_length,percent_identity,gap_openings,mismatches,e_value,bit_score,coverage"
    diamond blastp \
        --query "${queries}" \
        -d "${db}" \
        --ultra-sensitive \
        --evalue 1e-3 \
        --outfmt 6 qseqid sseqid qlen slen length pident gapopen mismatch evalue bitscore \
        --max-target-seqs 25 \
        --threads ${task.cpus} \
    | sort -k6 -n \
    | awk -v OFS='\\t' '
        BEGIN { print "'"\${header//,/\$'\\t'}"'" } 
        { split(\$1, qid, "|"); split(\$2, tid, "|"); print qid[2], tid[2], \$0, \$5/\$3 }
    ' \
    > hits.tsv

    # rough histogram
    awk -F'\\t' -v OFS='\\t' '
        BEGIN { print "species_target_uniprot_id","ortholog_count" }
        (NR == 1) { for ( i=0; i<=NF; i++ ) a[\$i]=i }
        (NR > 1) { c[\$a["species_target_uniprot_id"]]++ }
        END { for ( p in c ) print p, c[p] }
    ' hits.tsv \
    > hit_count.tsv

    """

}