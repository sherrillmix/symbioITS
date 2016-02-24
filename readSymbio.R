source("~/scripts/R/dna.R")
library(levenR)
library(dnaplotr)

if(!file.exists('SRR1257355.sra'))system('wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR125/SRR1257355/SRR1257355.sra')
if(!file.exists('SRR1257355.fastq.gz'))system('fastq-dump -gzip SRR1257355.sra')

if(!exists('x'))x<-read.fastq('SRR1257355.fastq.gz')
regex<-'^TCAG([ACGT]{8})GAATTGCAGAACTCCGTG' #CCATCTCATCCCTGCGTGTCTCCGAC__TCAG(N)8GAATTGCAGAACTCCGTG
selector<-grepl(regex,x$seq)
symbio<-x[selector,]
symbio$barcode<-sub(sprintf('%s.*',regex),'\\1',symbio$seq)
bigBarcodes<-table(symbio$barcode)
bigBarcodes<-bigBarcodes[bigBarcodes>3000]
if(length(bigBarcodes)!=15)stop(simpleError('Did not find 15 samples'))
symbio<-symbio[symbio$barcode %in% names(bigBarcodes),]
symbio$trim<-sub('N.*','',sub(regex,'',symbio$seq))
symbio<-symbio[nchar(symbio$trim)>300&nchar(symbio$trim<400),]
#http://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=244965
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4285332/table/tbl1/
barcodes<-read.csv('barcodes.csv',stringsAsFactors=FALSE)
if(any(sort(names(bigBarcodes))!=sort(barcodes$barcode)))stop(simpleError('Mismatch between observed and expected barcodes'))
rownames(barcodes)<-barcodes$barcode
symbio$sample<-barcodes[symbio$barcode,'sample']


bigSeqProp<-tapply(symbio$trim,symbio$barcode,function(x)max(table(x))/length(x))
bigSeqs<-tapply(symbio$trim,symbio$barcode,function(x)tail(names(sort(table(x))),1))
aligns<-cacheOperation('work/align.Rdat',lapply,barcodes$barcode,function(bar){
})

for(bar in unique(symbio$barcode)){
	png(sprintf('out/bar%s.png',bar),width=3000,height=1500,res=250)
		plotDNA(sort(symbio[symbio$barcode==bar,'trim']),main=with(barcodes[bar,],sprintf('%s (%s %s %s)',sample,host,dgge,ifelse(isoclonal,'isoclonal','field'))))
	dev.off()
}
