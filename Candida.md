
Candida glabrata 属于酵母科中酵母属，而其他三种念珠菌Candida tropicalis、Candida albicans、Candida parapsilosis属于德巴利酵母科念珠菌属


```shell

nwr member Candida |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

nwr lineage Candida |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

nwr info Nakaseomyces

# Nakaseomyces - genus
# --------------------
# NCBI Taxonomy ID: 374468
# First description:
# * Nakaseomyces Kurtzman, 2003
# Part of the Plants and Fungi.

# Comments: code compliant
```

| rank | count |
| --- | ---:|
| genus | 1 |
| species | 48 |
| no rank | 1 |
| strain | 111 |
| isolate | 1 |


| #rank | sci_name | tax_id |
| --- | --- | --- |
| kingdom | Fungi | 4751 |
| subkingdom | Dikarya | 451864 |
| phylum | Ascomycota | 4890 |
| subphylum | Saccharomycotina | 147537 |
| class | Saccharomycetes | 4891 |
| order | Saccharomycetales | 4892 |
| family | Debaryomycetaceae | 766764 |
| *genus* | Candida | 5475 |


```shell

mkdir -p ~/genome/Candida/summary
cd ~/genome/Candida/summary

nwr member Debaryomycetaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

echo -e "374468\tNakaseomyces\tgenus\tPlants and Fungi" >> genus.list.tsv

# 21 genus.list.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %' -- Crossbreeding of two species
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done


```
Download all assemblies

```shell
cd ~/genome/Aspergillus/summary

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Saccharomyces')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv


SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv

SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv

cat raw.tsv |
    tsv-uniq |
    datamash check

# 378 lines, 7 fields



cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-6 |
    perl ~/genome/Scripts/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[6]}++; # abbr_name
        $seen{$F[6]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\t%s\n}, $F[6], $F[3], $F[4], $F[1], $F[5];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- tsv-sort -k4,4 -k1,1 \
    > Candida.assembly.tsv

datamash check < Candida.assembly.tsv
# 377 lines, 5 fields
cat Candida.assembly.tsv | tsv-uniq -f 1 --repeated

cat Candida.assembly.tsv | tsv-filter --str-not-in-fld 2:ftp


rm raw*.*sv

```
Count before download

```shell

cd ~/genome/Candida

nwr template summary/Candida.assembly.tsv \
    --count \
    --rank genus


bash Count/strains.sh

bash Count/rank.sh

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| genus | #species | #strains |
|---|--:|--:|
| Aciculoconidium | 1 | 1 |
| Candida | 18 | 186 |
| Danielozyma | 1 | 1 |
| Debaryomyces | 7 | 21 |
| Hyphopichia | 6 | 8 |
| Kodamaea | 5 | 11 |
| Kurtzmaniella | 8 | 8 |
| Lodderomyces | 1 | 3 |
| Meyerozyma | 4 | 23 |
| Millerozyma | 2 | 3 |
| Nakaseomyces | 6 | 62 |
| Priceomyces | 4 | 4 |
| Saccharomyces | 1 | 1 |
| Scheffersomyces | 7 | 10 |
| Spathaspora | 9 | 11 |
| Suhomyces | 4 | 4 |
| Teunomyces | 3 | 3 |
| Wickerhamia | 1 | 1 |
| Yamadazyma | 12 | 15 |


Download and check

```shell

nwr template summary/Candida.assembly.tsv --ass

bash ASSEMBLY/rsync.sh

bash ASSEMBLY/check.sh

bash ASSEMBLY/n50.sh 100000 1000 1000000

cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"

# N50_min S_min   C_max
# 59311   10609954        536



cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
# S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
# 10902725.5      13049481        29999   500403  113.5   1409.5

bash ASSEMBLY/collect.sh

bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```

| #item | fields | lines |
|---|--:|--:|
| url.tsv | 3 | 376 |
| check.lst | 1 | 376 |
| collect.tsv | 20 | 377 |
| n50.tsv | 4 | 377 |
| n50.pass.tsv | 4 | 299 |
| collect.pass.tsv | 23 | 299 |
| pass.lst | 1 | 298 |
| omit.lst | 1 | 283 |
| rep.lst | 1 | 93 |

Rsync to hpcc



```shell



# <!-- 371 lines, 33 fields -->


nwr template summary/Candida.assembly.tsv \
    --mh \
    --parallel 48 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4







cat MinHash/abnormal.lst
# De_han_C11_GCA_016097515_1
# De_han_J6_GCA_001682995_1
# Me_cari_MG20W_GCA_000755205_1
# N_gla_15_0126687_GCA_026540015_1
# N_gla_16_0107051_GCA_026539545_1
# N_gla_16_0452821_GCA_026539665_1
# N_gla_16_6715916_GCA_026539885_1
# N_gla_16_6793165_GCA_026539985_1
# N_gla_17_6506958_GCA_026539955_1
# N_gla_17_6789296_GCA_026539265_1
# Sc_she_ATY839_GCA_002118035_1
# Sc_she_NBRC_1983_GCA_002118155_1



#!/usr/bin bash

ARRAY=('species::2')

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/genome/Scripts/condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.species.newick | rsvg-convert -o Trichoderma.minhash.png



```

| genus | #species | #strains |
|---|--:|--:|
| Aciculoconidium | 1 | 1 |
| Candida | 18 | 186 |
| Danielozyma | 1 | 1 |
| Debaryomyces | 7 | 21 |
| Hyphopichia | 6 | 8 |
| Kodamaea | 5 | 11 |
| Kurtzmaniella | 8 | 8 |
| Lodderomyces | 1 | 3 |
| Meyerozyma | 4 | 23 |
| Millerozyma | 2 | 3 |
| Nakaseomyces | 6 | 62 |
| Priceomyces | 4 | 4 |
| Saccharomyces | 1 | 1 |
| Scheffersomyces | 7 | 10 |
| Spathaspora | 9 | 11 |
| Suhomyces | 4 | 4 |
| Teunomyces | 3 | 3 |
| Wickerhamia | 1 | 1 |
| Yamadazyma | 12 | 15 |

| #family | genus | species | count |
| --- | --- | --- | ---:|
| Debaryomycetaceae | Candida | Candida albicans | 101 |
| Saccharomycetaceae | Nakaseomyces | Nakaseomyces glabratus | 56 |


| genus | #species | #strains |
|---|--:|--:|
| Candida | 14 | 48 |
| Debaryomyces | 1 | 1 |
| Hyphopichia | 1 | 1 |
| Kurtzmaniella | 1 | 1 |
| Lodderomyces | 1 | 1 |
| Meyerozyma | 1 | 1 |
| Millerozyma | 1 | 1 |
| Nakaseomyces | 1 | 25 |
| Saccharomyces | 1 | 1 |
| Scheffersomyces | 2 | 3 |
| Spathaspora | 2 | 2 |
| Suhomyces | 1 | 1 |
| Yamadazyma | 2 | 3 |

| #item | count |
| --- | --- |
| Proteins | 524,772 |
| Unique headers and annotations | 524,772 |
| Unique proteins | 524,772 |
| all.replace.fa | 524,772 |
| all.annotation.tsv | 524,773 |
| all.info.tsv | 524,773 |


ARRAY=(
    'Can_alb::C_alb_19F_GCA_000775445_1' # 1
    'Can_tro::C_tro_MYA_3404_GCA_013177555_1' # 5
    'Can_par::C_par_90_137_GCA_011316035_2' # 4
    'Can_gla::N_gla_040_PSC_GCA_024666015_1' # 19
)


| #Serial | Group | Target | Count |
| --- | --- | --- | ---:|
| 1 | Can_alb | C_alb_19F_GCA_000775445_1 | 60 |
| 5 | Can_tro | C_tro_MYA_3404_GCA_013177555_1 | 7 |
| 4 | Can_par | C_par_90_137_GCA_011316035_2 | 29 |
| 19 | Can_gla | N_gla_040_PSC_GCA_024666015_1 | 53 |





























