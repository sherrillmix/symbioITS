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
symbio<-symbio[nchar(symbio$trim)>300&nchar(symbio$trim)<400,]
#http://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=244965
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4285332/table/tbl1/
barcodes<-read.csv('barcodes.csv',stringsAsFactors=FALSE)
if(any(sort(names(bigBarcodes))!=sort(barcodes$barcode)))stop(simpleError('Mismatch between observed and expected barcodes'))
rownames(barcodes)<-barcodes$barcode
symbio$sample<-barcodes[symbio$barcode,'sample']


bigSeqProp<-tapply(symbio$trim,symbio$barcode,function(x)max(table(x))/length(x))
bigSeqs<-tapply(symbio$trim,symbio$barcode,function(x)tail(names(sort(table(x))),1))
aligns<-cacheOperation('work/align.Rdat',lapply,barcodes$barcode,function(bar,bigSeqs){
	message(bar)
	align<-levenAlign(symbio[symbio$barcode==bar,'trim'],bigSeqs[bar],nThreads=10,homoLimit=3,append=c(TRUE,TRUE))
	dist<-leven(symbio[symbio$barcode==bar,'trim'],bigSeqs[bar],nThreads=10,homoLimit=3,append=c(TRUE,TRUE))[,1]
	return(c(align,'dist'=list(dist)))
},bigSeqs)
names(aligns)<-barcodes$barcode

uniqSeq<-table(symbio$trim)
write.fa(paste(1:length(uniqSeq),uniqSeq,sep='_'),names(uniqSeq),'work/uniq.fa')
if(!file.exists('work/uniq.swarm'))system('~/installs/swarm/swarm work/uniq.fa -f -o work/uniq.swarm')

otus<-lapply(strsplit(readLines('work/uniq.swarm'),' '),function(x)as.numeric(sub('_[0-9]+$','',x)))
otuLookup<-data.frame('id'=unlist(otus),'otu'=rep(1:length(otus),sapply(otus,length)))
rownames(otuLookup)<-names(uniqSeq)[otuLookup$id]
symbio$otu<-otuLookup[symbio$trim,'otu']

otuN<-tapply(symbio$trim,symbio$otu,function(x)length(x))
otuSeq<-tapply(symbio$trim,symbio$otu,function(x)tail(names(sort(table(x))),1))

otuTable<-table(symbio$sample,symbio$otu)
bigOtus<-otuTable[,names(otuN)[otuN>15]]
bigOtuDist<-leven(otuSeq[colnames(bigOtus)],append=2,homoLimit=3,nThreads=10)
bigOtuProps<-t(apply(bigOtus,1,function(x)x/sum(x)))
rownames(bigOtuProps)<-sapply(rownames(bigOtuProps),function(x)with(barcodes[barcodes$sample==x,],sub('CCMP2467:rt-147 ','B1:A1',paste(sub('([A-Z])[a-z]+ ','\\1. ',host),dgge,ifelse(isoclonal,'clone','field')))))
clust<-hclust(as.dist(bigOtuDist))
cols<-c('white',rainbow.lab(99))
breaks<-c(-1e-9,10^seq(log10(min(bigOtuProps[bigOtuProps>0])-1e-9),log10(1+1e-9),length.out=100))
pdf('out/tree.pdf')
	heatmap(bigOtuProps,col=cols,breaks=breaks,scale='none',Colv=as.dendrogram(clust),mar=c(4,11),labCol='',xlab='',add.expr={box();abline(v=.5+1:100,h=.5+1:100,col='#00000022');mtext('OTU cluster',1,.4)})
	prettyLabs<-unique(round(pretty(range(log10(breaks[-1])))))
	prettyLabs<-prettyLabs[10^prettyLabs>min(breaks[-1])]
	insetPos<-c(grconvertX(0.015,'nfc','user'),grconvertY(0.025,'nfc','user'),grconvertX(0.4,'nfc','user'),grconvertY(0.04,'nfc','user'))
  breakPos<-(log10(breaks[-1])-log10(min(breaks[-1])))/max(log10(breaks[-1])-log10(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-1]+1e-3,insetPos[2],breakPos[-length(breakPos)],insetPos[4],col=cols[-1],xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyPos<-(prettyLabs-log10(min(breaks[-1])))/(log10(max(breaks[-1]))-log10(min(breaks[-1])))*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,sub('\\.?0+$','',format(10^prettyLabs,scientific=FALSE)),xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.45,"Proportion of sample",xpd=NA,adj=c(.5,0))
dev.off()



#distMat<-cacheOperation('work/bigMat.Rdat',leven,uniqSeq,nThreads=12,homoLimit=3,append=2,vocal=100)

for(ii in 1:nrow(barcodes)){
	bar<-barcodes[ii,'barcode']
	message(bar)
	png(sprintf('out/strain%02d_%s.png',ii,barcodes[bar,'sample']),width=3000,height=1500,res=250)
		#plotDNA(sort(symbio[symbio$barcode==bar,'trim']),main=with(barcodes[bar,],sprintf('%s (%s %s %s)',sample,host,dgge,ifelse(isoclonal,'isoclonal','field'))))
		#seqs<-aligns[[bar]][[2]][order(aligns[[bar]][[3]],aligns[[bar]][[2]])]
		seqs<-aligns[[bar]][[2]]
		groups<-symbio[symbio$barcode==bar,'otu']
		dists<-aligns[[bar]][[3]]
		groupDists<-tapply(dists,groups,mean)
		groupRanks<-rank(groupDists,ties.method='first')
		plotDNA(removeGapCols(seqs[order(groupDists[as.character(groups)],groups,seqs)],maxGapProp=.95),main=with(barcodes[bar,],sprintf('%s (%s %s %s)',sample,host,dgge,ifelse(isoclonal,'isoclonal','field'))))
	dev.off()
}

cladeOtus<-as.numeric(tail(names(sort(apply(otuTable,2,sum))),3))
otuAligns<-cacheOperation('work/otuAlign.Rdat',lapply,cladeOtus,function(otu,bigSeqs){
	message(otu)
	thisSeqs<-symbio[symbio$otu==otu,'trim']
	otuSeq<-tail(names(sort(table(thisSeqs))),1)
	align<-levenAlign(thisSeqs,otuSeq,nThreads=12,homoLimit=3,append=c(TRUE,TRUE))
	dist<-leven(thisSeqs,otuSeq,nThreads=12,homoLimit=3,append=c(TRUE,TRUE))[,1]
	return(c(align,'dist'=list(dist)))
},bigSeqs)
names(otuAligns)<-cladeOtus

for(ii in cladeOtus){
	message(ii)
	png(sprintf('out/otu%02ds.png',ii),width=3000,height=1500,res=250)
		seqs<-otuAligns[[as.character(ii)]][[2]]
		groups<-symbio[symbio$otu==ii,'sample']
		dists<-otuAligns[[as.character(ii)]][[3]]
		groupDists<-tapply(dists,groups,mean)
		clade<-tail(names(sort(table(symbio[symbio$otu==ii,'dgge']))),1)
		groupRanks<-rank(groupDists,ties.method='first')
		plotDNA(removeGapCols(seqs[order(groupDists[as.character(groups)],groups,seqs)],maxGapProp=.95),main=sprintf('OTU %d (clade %s)',ii,clade),groups=groups)
	dev.off()
}

