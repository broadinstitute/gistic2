# 2017-03-21  obtain and parse GENCODE v22 gene files for GISTIC refgene

# get gencode gene annotation file
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
gunzip gencode.v22.annotation.gtf.gz
# extract genes as table (with attributes as columns)
cat gencode.v22.annotation.gtf | perl parse_gencode_genes.pl > gencode_genes.tsv

# get Entrez-to-ENSEMBLE identifier mapping from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Homo_sapiens.gene_info.gz
# extract table
cat Homo_sapiens.gene_info | perl parse_ncbi2ensembl.pl > ncbi2ensembl.tsv


