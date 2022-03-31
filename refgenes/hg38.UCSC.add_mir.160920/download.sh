wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz -O cytoBand.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz -O refGene.txt.gz
## NOTE: download refLink.txt and refSeqStatus.txt from the table browser - can't find it in UCSC ftp or http downlaod areas
##wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refLink.txt.gz -O refLink.txt.gz
##wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refSeqStatus.txt.gz -O refSeqStatus.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refLink.txt.gz -O refLink.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgRna.txt.gz -O wgRna.txt.gz
gunzip -N *.txt.gz
