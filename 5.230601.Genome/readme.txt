*.dmp files are bcp-like dump from GenBank taxonomy database.

General information.
Field terminator is "\t|\t"
Row terminator is "\t|\n"

nodes.dmp file consists of taxonomy nodes. The description for each node includes the following
fields:
        tax_id                                  -- node id in GenBank taxonomy database
        parent tax_id                           -- parent node id in GenBank taxonomy database
        rank                                    -- rank of this node (superkingdom, kingdom, ...)
        embl code                               -- locus-name prefix; not unique
        division id                             -- see division.dmp file
        inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
        genetic code id                         -- see gencode.dmp file
        inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
        mitochondrial genetic code id           -- see gencode.dmp file
        inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
        GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
        hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
        comments                                -- free-text comments and citations

Taxonomy names file (names.dmp):
        tax_id                                  -- the id of node associated with this name
        name_txt                                -- name itself
        unique name                             -- the unique variant of this name if name not unique
        name class                              -- (synonym, common name, ...)

Divisions file (division.dmp):
        division id                             -- taxonomy database division id
        division cde                            -- GenBank division code (three characters, e.g. BCT, PLN, VRT, MAM, etc.)
        division name                           -- e.g. Bacteria, Plants and Fungi, Vertebrates, Mammals, etc.
        comments

Genetic codes file (gencode.dmp):
        genetic code id                         -- GenBank genetic code id
        abbreviation                            -- genetic code name abbreviation
        name                                    -- genetic code name
        cde                                     -- translation table for this genetic code
        starts                                  -- start codons for this genetic code

Deleted nodes file (delnodes.dmp):
        tax_id                                  -- deleted node id

Merged nodes file (merged.dmp):
        old_tax_id                              -- id of nodes which has been merged
        new_tax_id                              -- id of nodes which is result of merging

Citations file (citations.dmp):
        cit_id                                  -- the unique id of citation
        cit_key                                 -- citation key
        pubmed_id                               -- unique id in PubMed database (0 if not in PubMed)
        medline_id                              -- unique id in MedLine database (0 if not in MedLine)
        url                                     -- URL associated with citation
        text                                    -- any text (usually article name and authors).
                                                -- The following characters are escaped in this text by a backslash:                                                -- newline (appear as "\n"),
                                                -- tab character ("\t"),
                                                -- double quotes ('\"'),
                                                -- backslash character ("\\").
        taxid_list                              -- list of node ids separated by a single space

Organism images file (images.dmp):
        image_id                                -- the unique id of image
        image_key                               -- image key
        url                                     -- image URL associated with citation
        license                                 -- image license
        attribution                             -- image attribution
        source                                  -- source of the image
        properties                              -- various image properties separated by semicolon
        taxid_list                              -- list of node ids separated by a single space









*.dmp文件是来自GenBank分类数据库的类似bcp的转储。

一般信息。
字段终止符是"\t|\t"
行终止符是"\t|n"

nodes.dmp文件由分类学节点组成。每个节点的描述包括以下内容
字段：
        tax_id -- GenBank分类学数据库中的节点ID
        parent tax_id -- GenBank分类学数据库中的父节点ID
        rank -- 这个节点的等级（superkingdom, kingdom, ...）。
        embl code -- 定位名称的前缀；不是唯一的
        division id -- 参见division.dmp文件
        inherited div flag (1 or 0) -- 如果节点从父本继承了分部，则为1
        遗传代码ID -- 见gencode.dmp文件
        继承的GC标志（1或0） -- 如果节点从父母那里继承了遗传代码，则为1
        线粒体遗传密码ID -- 见gencode.dmp文件
        继承的MGC标志（1或0） -- 如果节点从父母那里继承了线粒体基因码，则为1
        GenBank隐藏标志（1或0） -- 如果名字在GenBank条目线中被抑制，则为1
        Hidden subree root flag (1 or 0) -- 如果这个子树还没有序列数据，则为1
        评论 -- 自由文本的评论和引用

分类学名称文件（names.dmp）：
        tax_id -- 与该名称相关的节点的id
        name_txt -- 名称本身
        唯一的名字 -- 如果名字不是唯一的，则为该名字的唯一变体
        名称类别 -- （同义词，通用名称，...）。
分区文件（division.dmp）：
        Division id -- 分类数据库的分部ID
        division cde -- GenBank分部代码（三个字符，如BCT、PLN、VRT、MAM等）。
        分部名称 -- 例如：细菌、植物和真菌、脊椎动物、哺乳动物等。
        评论

遗传密码文件（gencode.dmp）：
        遗传密码ID -- GenBank的遗传密码ID
        缩写 -- 遗传密码名称缩写
        名称 -- 遗传密码名称
        cde -- 该遗传密码的翻译表
        starts -- 该遗传密码的起始密码子

删除的节点文件（delnodes.dmp）：
        tax_id -- 删除的节点ID

合并的节点文件（merged.dmp）：
        old_tax_id -- 已经合并的节点的id
        new_tax_id -- 合并后的节点的id。

引用文件（citations.dmp）：
        cit_id -- 引文的唯一ID
        cit_key -- 引文关键词
        pubmed_id -- 在PubMed数据库中的唯一ID（如果不在PubMed中则为0）
        medline_id -- MedLine数据库中的唯一ID（如果不在MedLine中则为0）
        url -- 与引文相关的URL
        text -- 任何文本（通常是文章名称和作者）。
                                                -- 以下字符在此文本中用反斜杠转义：
                                                --换行符（显示为"\n"）、
                                                --制表符（"\t"）、
                                                --双引号（"\"）、
                                                -- 反斜杠字符（"\"）。
        taxid_list -- 节点ID的列表，用一个空格分隔

生物体图像文件（images.dmp）：
        image_id -- 图像的唯一ID
        image_key -- 图像密钥
        url -- 与引文相关的图像URL
        license -- 图像许可证
        attribution -- 图像的归属
        source -- 图像的来源
        properties -- 各种图像属性，用分号分隔
        taxid_list -- 节点ID的列表，用一个空格分隔