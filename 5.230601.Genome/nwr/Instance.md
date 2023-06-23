nwr的一个使用流程

### 准备工作

```shell
brew install wang-q/tap/nwr
brew install wang-q/tap/tsv-utils
brew install sqlite
brew install miller


nwr download

nwr txdb

nwr ardb
nwr ardb --genbank


```

### NCBI ASSEMBLY

* assembly_level

```shell

for i in refseq genbank; do
    cat ~/.nwr/assembly_summary_${C}.txt | sed '1d' | 
        tsv-summarize -H -g assembly_level,genome_rep --count |
        keep-header -- sort |
        mlr --itsv --omd cat |   # 这一步是让输出格式符合markdown的三线表格式
        sed 's/-\s*|$/-:|/'| echo -e "\nTable: ${C}\n\n"
done

# 统计NCBI上所有物种的装配水平的数量

```
| assembly_level | genome_rep | count |
| --- | --- | ---:|
| Chromosome | Full | 5635 |
| Chromosome | Partial | 915 |
| Complete Genome | Full | 47937 |
| Complete Genome | Partial | 25 |
| Contig | Full | 157432 |
| Contig | Partial | 1 |
| Scaffold | Full | 93586 |

Table: refseq

#### Example 1: count qualified assemblies of Eukaryote groups（以Ascomycetes::Ascomycota为例）

```shell

# 查看Ascomycota的属水平的的装配情况



echo -e "GROUP_NAME\tSCI_NAME\tComplete Genome\tChromosome\tScaffold\tContig" \
    > groups.tsv

GROUP_NAME=Ascomycetes
SCI_NAME=Ascomycota



GENUS=$(
    nwr member Ascomycota -r genus |
        grep -v -i "Candidatus " | #排除未命名的属
        grep -v -i "candidate " |
        sed '1d' | #去除表头，表头内容为：#tax_id   sci_name    rank    division
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$/\)/' |
        sed 's/^/\(/'    # 将输出的列表的第一行储存在一个数组中
)


printf "$GROUP_NAME\t$SCI_NAME\t"

for L in 'Complete Genome' 'Chromosome' 'Scaffold' 'Contig'; do
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND genus_id IN $GENUS
                AND assembly_level IN ('$L')
            " |   # 这个地方应该是SQLite数据库的格式
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    done |
    tr "\n" "\t" |
    sed 's/\t$//'

    echo;
done \
    >> groups.tsv

cat groups.tsv |
    mlr --itsv --omd cat

```

#### Example 2：find accessions of a species（Staphylococcus capitis - 29388 - 头状葡萄球菌）

下面的步骤中涉及到SQLite的使用，可以参考[此页面](./SQLite.md)

```shell
nwr info "Staphylococcus capitis"

# taphylococcus capitis - species
# --------------------------------
# NCBI Taxonomy ID: 29388
# First description:
# * Staphylococcus capitis Kloos and Schleifer 1975
# Part of the Bacteria.

nwr member 29388

# tax_id sci_name        rank    division
# 29388   Staphylococcus capitis  species Bacteria
# 1441378 Staphylococcus capitis AYP1020  strain  Bacteria
# 1296619 Staphylococcus capitis CR01     strain  Bacteria
# 1189311 Staphylococcus capitis QN1      strain  Bacteria
# 904334  Staphylococcus capitis VCU116   strain  Bacteria
# 553212  Staphylococcus capitis SK14     strain  Bacteria
# 435838  Staphylococcus capitis C87      strain  Bacteria
# 74703   Staphylococcus capitis subsp. urealyticus       subspecies      Bacteria
# 72758   Staphylococcus capitis subsp. capitis   subspecies      Bacteria

```

```shell
echo '
.headers ON

SELECT
    organism_name,
    species,
    genus,
    ftp_path,
    assembly_level
FROM ar
WHERE 1=1
    AND tax_id != species_id    -- with strain ID
    AND species_id IN (29388)
' |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > Scap.assembly.tsv

echo '
SELECT
    species || " " || REPLACE(assembly_accession, ".", "_") AS organism_name,
    species,
    genus,
    ftp_path,
    assembly_level
FROM ar
WHERE 1=1
    AND tax_id = species_id     -- no strain ID
    AND assembly_level IN ("Chromosome", "Complete Genome")
    AND species_id IN (29388)
' |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> Scap.assembly.tsv



```
#### Example 3：find model organisms in a family

```shell
echo "
.headers ON

    SELECT
        tax_id,
        organism_name
    FROM ar
    WHERE 1=1
        AND family IN ('Enterobacteriaceae')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    sed '1s/^/#/' |
    mlr --itsv --omd cat



```





cat mash.dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.5) # k=3
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '





