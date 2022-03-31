require(data.table)
require(GenomicRanges)

RG <- fread("gencode_genes.tsv")
nrow(RG) #=> 60483

sort(table(RG$gene_type),decreasing=TRUE)
##=> protein_coding: 19814
##   processed_pseudogene: 10304
##   lincRNA: 7656
##   antisense: 5565
##   miRNA: 4093
##   unprocessed_pseudogene: 2574
##   misc_RNA: 2298
##   ...


dupnames <- RG$gene_name[duplicated(RG$gene_name)]
length(dupnames) #=> 2096
length(unique(dupnames)) #=> 399

RGD <- RG[RG$gene_name %in% dupnames,]
nrow(RGD) #=> 2495

sort(table(RGD$gene_name),decreasing=TRUE)
##=> Y_RNA: 782
##   Metazoa_SRP: 306
##   U3: 56
##   uc_338: 36
##   U6: 34
##   snoU13: 34
##   ...

table(RGD$gene_type[grepl("^SNORD[0-9]+",RGD$gene_name)])
##=> snoRNA: 142
table(RGD$gene_type[grepl("^SNORA[0-9]+",RGD$gene_name)])
##=> snoRNA: 364
sort(table(RGD$gene_type[RGD$gene_name=="Y_RNA"]))
##=> misc_RNA: 782
##=> snoRNA: 653
##=> protein_coding: 241
##=> snRNA: 102

## count duplications of each name
RGD[,num.dup := .N,by=gene_name]

## count overlaps
RGD.gr <- GRanges(seqnames=RGD$chr,ranges=IRanges(start=RGD$start,end=RGD$end),strand=RGD$strand,gene_id=RGD$gene_id,gene_name=RGD$gene_name)

lap.count <- function(gene_name) {
    gr <- RGD.gr[RGD.gr$gene_name==gene_name]
    olaps <- findOverlaps(gr,gr)
    max(as.vector(table(queryHits(olaps))))
}

RGD[,num.lap := lap.count(gene_name),by=gene_name]


xy_dupnames <- RGD$gene_name[grep("ENSGR",RGD$gene_id)]
RGD[,xy := as.numeric(gene_name %in% xy_dupnames & chr %in% c("chrX","chrY")) ]


write.tsv <- function(x,file,...) { write.table(x,file,sep="\t",quote=FALSE,row.names=FALSE,...) }
write.tsv(RGD,file="duplicate_gencode_names.170322.txt")

write.table(table(RG$gene_type,RG$level),"type_vs_level.170322.tsv",sep="\t",quote=FALSE,row.names=TRUE)

UCSC <- fread("/xchip/gistic/variables/hg38/hg38.UCSC.add_miR.160920.refgene.txt")
RG[,in.UCSC := RG$gene_name %in% UCSC$symb]

write.table(table(RG$gene_type[RG$in.UCSC],RG$level[RG$in.UCSC]),"type_vs_level.in_UCSC.170322.tsv",sep="\t",quote=FALSE,row.names=TRUE)

SCNA <- fread("/xchip/gistic/schum/ICGC/nuria_161006/driver_genes/scna_drivers.161206.txt")
RG[,scna.driver := RG$gene_name %in% SCNA$gene]

