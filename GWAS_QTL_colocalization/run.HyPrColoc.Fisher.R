df=data.frame(t(rep(NA,5)))
m=read.table('results/eQTL.source.counts.0.98.txt');
n=read.table(pipe('grep -P \'\tge\t\' data/external_eQTLs/tabix_ftp_paths.tsv | grep -v GTEx | awk -F \'\t\' \'{print $1"_ge_"$2"\t"$8}\''),row.names=1)
n=rbind(n,1367);rownames(n)[nrow(n)]="i2QTL"
n=rbind(n,10000);rownames(n)[nrow(n)]="eQTLGen"
n$V2=n[as.character(m$V2),]
rownames(n)=as.character(m$V2)
ss=tapply(n$V2,m$V1,mean)

for (typ in unique(m$V1)) {
	print(typ);
	ft=fisher.test(cbind(t(cbind(sum(subset(m,V1%in%typ,'V3')),sum(m$V3)-sum(subset(m,V1%in%typ,'V3')))),t(cbind(sum(subset(m,V1%in%typ,'V4')),sum(m$V4)-sum(subset(m,V1%in%typ,'V4'))))));ft2=fisher.test(cbind(t(cbind(sum(subset(m,V1%in%typ,'V3')),sum(m$V3)-sum(subset(m,V1%in%typ,'V3')))),t(cbind(sum(subset(m,V1%in%typ,'V4')),sum(m$V4)-sum(subset(m,V1%in%typ,'V4'))))),alternative='greater');df=rbind(df,c(typ,ft$estimate,ft$conf.int[1],ft$conf.int[2],ft$p.value))
}
df=df[-1,]
df=cbind(df,p.adjust(df$X5,method='bonferroni'))
df=df[order(df$X2,decreasing=T),]
df
cor.test(as.numeric(df$X2),ss[df$X1])
