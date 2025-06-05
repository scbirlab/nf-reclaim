process fetch_gnomad_constraints {

    tag "v${version}"

    publishDir( 
        "${params.outputs}/gnomad", 
        mode: 'copy',
        saveAs: { "${it}" }
    )

    input:
    val version

    output:
    path "*.constraint_metrics_.tsv"

    script:
    """
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/${version}/constraint/gnomad.v${version}.constraint_metrics.tsv

    python -c '
    import pandas as pd
    transcript = "protein_coding"
    (   
        pd.read_csv("gnomad.v${version}.constraint_metrics.tsv", sep="\t")
        .query("gene.str.len() > 0")
        .query("canonical and mane_select and transcript_type == @transcript")
        .assign(taxon_id=9606)
        [["taxon_id", "gene", "lof_hc_lc.pLI", "lof.oe_ci.upper"]]
        .rename(columns={"lof_hc_lc.pLI": "pLI", "lof.oe_ci.upper": "LOEUF"})
        .groupby(["taxon_id", "gene"])
        .apply(lambda x: x.nsmallest(1, "LOEUF"), include_groups=False)
        .reset_index()
        .to_csv("gnomad.v${version}.constraint_metrics_.tsv", sep="\t", index=False)
    )
    '

    """

}