process fetch_pubchem_vendors {

    tag "${id}"

    // errorStrategy 'retry'
    // maxRetries 2

    publishDir( 
        "${params.outputs}/inhibitors", 
        mode: 'copy',
        saveAs: { "${id}-${it}" }
    )

    input:
    tuple val( id ), path( table )
    val pubchem_url

    output:
    tuple val( id ), path( "vendors.tsv" )

    script:
    """
    set -x
    echo "=== DNS ==="
    getent hosts pubchem.ncbi.nlm.nih.gov
    echo "\\n=== ENV PROXY ==="
    env | grep -i proxy || echo "no proxy vars set"

    id_col=\$(head -n1 "${table}" | tr \$'\\t' \$'\\n' | grep -n -Fw pubchem_id | cut -d: -f1)

    # look up vendors
    printf 'pubchem_id\\tvendor_name\\tcatalog_number\\tvendor_url\\n' > vendors.tsv
    for cid in \$(tail -n+2 "${table}" | cut -f"\$id_col")
    do
        curl -vk4 "${pubchem_url}/rest/pug_view/categories/compound/\$cid/JSON" > result.json
        jq -r '.SourceCategories.Categories[] | select( .Category == "Chemical Vendors" ).Sources[] | [.SourceName, .RegistryID, .SourceRecordURL] | @tsv' \
        < result.json \
        > info.tsv
        awk -v OFS='\\t' -v cid="\$cid" '{ print cid,\$0 }' info.tsv >> vendors.tsv
        sleep 0.1
    done

    head -n1 vendors.tsv | cat - <(tail -n+2 vendors.tsv | sort -k1 -n) > vendors-sorted.tsv
    head -n1 "${table}" | cat - <(tail -n+2 "${table}" | sort -k"\$id_col" -n) > pubchem_ids-sorted.tsv

    join --header -t\$'\\t' \
        -1 "\$id_col" -2 1 \
        pubchem_ids-sorted.tsv vendors-sorted.tsv \
    > vendors.tsv

    """

}