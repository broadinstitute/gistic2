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

write.table(RGD,file="duplicate_gencode_names.170321.txt",sep="\t",quote=FALSE,row.names=FALSE)
