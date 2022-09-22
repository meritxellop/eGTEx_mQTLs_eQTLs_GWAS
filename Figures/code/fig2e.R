# https://stackoverflow.com/questions/35175606/how-to-explode-donut-chart-slices-using-r-ggplot2
library(ggplot2)

data=read.table('data/counts.QTLcoloc.classes.txt',header=T,sep='\t')

#category	count	fraction	ymax	ymin	labelPosition	label
#A	19213	0.0671	0.0671	0.0	0.05	No eQTL overlap\n 19,213 (6.71%)
#B	61380	0.2145	0.2816	0.0671	0.14	Colocalized e/mQTL\n 61,380 (21.45%)
#C	108644	0.3797	0.6613	0.2816	0.40	Independent e/mQTL\n 108,644 (37.97%)
#D	96915	0.3387	1.0	0.6613	0.85	Uncertain colocalization\n 96,915 (33.87%)

# Create data.
dat <- data.frame(
  category=c("No eQTL overlap", "Colocalized e/mQTL", "Independent e/mQTL", "Uncertain colocalization"),
  count=c(19213, 61380, 108644, 96915)
)

dat$fraction = dat$count / sum(dat$count)
dat = dat[order(dat$fraction), ]
dat$ymax = cumsum(dat$fraction)-0.01
dat$ymin = c(0, head(dat$ymax, n=-1))+0.01
dat$all_=length(unique(dat$category))
dat$x1=dat$all_-(1:nrow(dat))*0.5+1
dat$x2=dat$all_-(1:nrow(dat))*0.5+2

p2=ggplot()+aes(ymin=0)+geom_rect(data=dat,aes(fill=category,ymax=ymax, ymin=ymin,  xmax=x1, xmin=x2),color = "white")+
  ylim(0,1)+
  xlim(c(0,3+length(unique(dat$category))))+
  coord_polar(theta="y", direction = -1) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(panel.border = element_blank())
p2
