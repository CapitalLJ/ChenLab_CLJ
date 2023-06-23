


```shell

nwr member Aspergillus |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

nwr lineage Aspergillus |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tFungi/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat


```

| rank | count |
| --- | ---:|
| genus | 1 |
| species | 612 |
| subgenus | 6 |
| strain | 92 |
| no rank | 2 |
| varietas | 37 |
| subspecies | 1 |
| isolate | 2 |

| #rank | sci_name | tax_id |
| --- | --- | --- |
| kingdom | Fungi | 4751 |
| subkingdom | Dikarya | 451864 |
| phylum | Ascomycota | 4890 |
| subphylum | Pezizomycotina | 147538 |
| class | Eurotiomycetes | 147545 |
| subclass | Eurotiomycetidae | 451871 |
| order | Eurotiales | 5042 |
| family | Aspergillaceae | 1131492 |
| *genus* | Aspergillus | 5052 |











```shell
mkdir -p ~/genome/Aspergillus/summary
cd ~/genome/Aspergillus/summary

nwr member Aspergillaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv


wc -l genus.list.tsv


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
    > Aspergillus.assembly.tsv


datamash check < Aspergillus.assembly.tsv

cat Aspergillus.assembly.tsv |
    tsv-uniq -f 1 --repeated

cat Aspergillus.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp


rm raw*.*sv
```



Count before download

```shell
cd ~/genome/Aspergillus

nwr template summary/Aspergillus.assembly.tsv \
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
| Aspergillus | 133 | 1072 |
| Evansstolkia | 1 | 1 |
| Hamigera | 1 | 1 |
| Monascus | 3 | 18 |
| Penicilliopsis | 1 | 1 |
| Penicillium | 103 | 459 |
| Saccharomyces | 1 | 1 |
| Xeromyces | 1 | 1 |



Download and check

```shell

nwr template summary/Aspergillus.assembly.tsv --asss

bash ASSEMBLY/rsync.sh

bash ASSEMBLY/check.sh

bash ASSEMBLY/n50.sh 100000 1000 1000000

cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" |
    tsv-summarize -H --min "N50,S" --max "C"

# N50_min S_min   C_max
# 66657   26053171        1115


cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"

# S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
# 28019617.5      32502288        95695.5 493788  341     1163.7


cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

```
| #item | fields | lines |
|---|--:|--:|
| url.tsv | 3 | 1,554 |
| check.lst | 1 | 1,554 |
| collect.tsv | 20 | 1,555 |
| n50.tsv | 4 | 1,555 |
| n50.pass.tsv | 4 | 1,287 |
| collect.pass.tsv | 23 | 1,287 |
| pass.lst | 1 | 1,286 |
| omit.lst | 1 | 869 |
| rep.lst | 1 | 217 |


BioSample

```shell

ulimit -n `ulimit -Hn`

nwr template summary/Aspergillus.assembly.tsv --bs

bash BioSample/download.sh

bash BioSample/collect.sh

datamash check < BioSample/biosample.tsv

# 1552 lines, 46 fields

```

Rsync to hpcc

```shell
rsync -rvP /home/llj/genome/Aspergillus/ wangq@202.119.37.251:llj/Aspergillus

```   

MinHash

```shell
cd ~/llj/Aspergillus

nwr template summary/Aspergillus.assembly.tsv \
    --mh \
    --parallel 16 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst
#T_har_CGMCC_20739_GCA_019097725_1
#T_har_Tr1_GCA_002894145_1
#T_har_ZL_811_GCA_021186515_1

# Non-redundant strains within species
bash MinHash/nr.sh

# Distances between all selected sketches, then hierarchical clustering
bash MinHash/dist.sh



```

Condense branches in the minhash tree

```shell
ARRAY=(
    'species::2'
)

rm minhash.condensed.map

CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/llj/Scripts/condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.species.newick |
    rsvg-convert -o Aspergillus.minhash.png



```

Count valid species and strains


For genomic alignments

```shell
cd ~/llj/Aspergillus

nwr template summary/Aspergillus.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh 50

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv



```


| genus | #species | #strains |
|---|--:|--:|
| Aspergillus | 110 | 839 |
| Hamigera | 1 | 1 |
| Monascus | 3 | 17 |
| Penicilliopsis | 1 | 1 |
| Penicillium | 94 | 386 |
| Saccharomyces | 1 | 1 |
| Xeromyces | 1 | 1 |


| #family | genus | species | count |
| --- | --- | --- | ---:|
| Aspergillaceae | Aspergillus | Aspergillus flavus | 188 |
|  |  | Aspergillus fumigatus | 212 |
|  |  | Aspergillus niger | 106 |
|  |  | Aspergillus oryzae | 102 |
|  | Penicillium | Penicillium chrysogenum | 78 |



For protein families

```shell

cd ~/llj/Aspergillus 

nwr template summary/Aspergillus.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus
```


| genus | #species | #strains |
|---|--:|--:|
| Aspergillus | 76 | 308 |
| Monascus | 1 | 2 |
| Penicilliopsis | 1 | 1 |
| Penicillium | 70 | 190 |
| Saccharomyces | 1 | 1 |


Collect proteins

```shell
cd ~/llj/Aspergillus

nwr template summary/Aspergillus.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat



```
| #item | count |
| --- | --- |
| Proteins | 5,335,421 |
| Unique headers and annotations | 5,335,421 |
| Unique proteins | 5,335,421 |
| all.replace.fa | 5,335,421 |
| all.annotation.tsv | 5,335,422 |
| all.info.tsv | 5,335,422 |

Phylogenetics with fungi61

Find corresponding proteins by hmmsearch

```shell



faops size Protein/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2

# 43750
# 24980




```

Groups and targets
```shell
ARRAY=(
    'A_acu::A_aculeati_GCA_016748025_1' # 1
    'A_ter::A_ter_45A_GCA_001630395_1' # 2
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_tam_pse::A_cae_CBS_763_97_GCF_009193585_1' # 4
    'A_amo_syd_vers::A_puu_AS33_GCA_026319425_1' # 5
    'A_chev_cri::A_cri_AC_1_GCA_020466005_1' # 6
    'A_flavus_ory_paras_soj::A_ory_TK_48_GCA_009687065_1' # 7
    'A_asp::A_asp_NRRL_4770_GCA_027569075_1' # 8
    'A_nig_wel::A_nig_CBS_113_50_GCA_023134475_1' # 9  #nig在group9和10中都存在
    'A_luc_nig_tub::A_luc_TK_87_GCA_016860325_1' # 10
    'A_cal::A_cal_FKII_L2_CM_DR1_GCA_022814525_1' # 11




    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'A_wes::A_mel_CBS_546_65_GCF_016097325_1' # 3
    'T_atrov_koningio::T_atrov_IMI_206040_GCF_000171015_1' # 6
    'T_cit_lon_ree::T_ree_QM6a_GCF_000167675_1' # 7
    'T_vire::T_vire_Gv29_8_GCF_000170995_1' # 8
)
```