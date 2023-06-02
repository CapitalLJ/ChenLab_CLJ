###  安装nwr并启动本地数据库

```shell
brew install wang-q/tap/nwr # 0.5.4 or above
brew install sqlite 
# SQLite是一种开源的嵌入式关系型数据库管理系统（DBMS），广泛应用于各种应用程序和平台中。

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```
关于nwr的说明，详见[此页面](https://github.com/wang-q/nwr)


# Bacteria
All genomes of Bacteria and Archaea, species by species.Download all genomes and analyze representative strains.

### List all ranks

域（Domain）、界（Kingdom）、门（Phylum）、纲（Class）、目（Order）、科（Family）、属（Genus）、种（Species）

```shell
nwr member Bacteria |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat


# 拉丁名中"sp."是某属中一未知种的意思，"spp."是某属中多个未知种的意思。 "sp."和"spp."是种species的缩写，"sp."是种的单数，"spp."是种的复数。
# tsv-summarize对标签分隔值文件中的字段进行聚合操作。

nwr member Archaea |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat



```