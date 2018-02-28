library(fastqcr)
library(ggplot2)
library(reshape2)

qc.dir <- "./"
qc <- qc_aggregate(qc.dir)
summary(qc)

ggplot(data.frame(x=as.numeric(qc$tot.seq)/1e6),aes(x=x))+theme_classic()+geom_histogram(binwidth=5)+scale_x_continuous('Total reads (M)')+ggtitle('Total reads / sample')+theme(plot.title=element_text(hjust=0.5))
ggsave('ReadsPerSample',height=4,width=4)

ggplot(data.frame(x=as.numeric(qc$pct.dup)),aes(x=x))+theme_classic()+geom_histogram(binwidth=5)+scale_x_continuous('% dup ')+ggtitle('Percent duplication / sample')+theme(plot.title=element_text(hjust=0.5))
ggsave('PctDupSample',height=4,width=4)

samples<-list.files(path = "./", pattern = "*.zip")
pbsq<-matrix(0,55,length(samples))
pbsc_g<-matrix(0,55,length(samples))
for (ii in 1: length(samples)) {
    qcf <- qc_read(samples[ii])
    pbsq[,ii]<-qcf$per_base_sequence_quality$Mean
    pbsc_g[,ii]<-qcf$per_base_sequence_content$G
}


ggplot(melt(pbsc_g),aes(x=Var1,y=value,colour=Var2,group=Var2))+theme_classic()+geom_line()+scale_y_continuous('Per base sequence content - G') +scale_x_continuous('')+theme(axis.text.x=element_blank(),legend.position='none')
ggsave('PerBaseSeqC_G1.pdf',height=4,width=4)

ggplot(melt(pbsq),aes(x=Var1,y=value,colour=Var2,group=Var2))+theme_classic()+geom_line()+scale_y_continuous('Per base sequence quality') +scale_x_continuous('')+theme(axis.text.x=element_blank(),legend.position='none')
ggsave('PerBaseSeqQ.pdf',height=4,width=4)
