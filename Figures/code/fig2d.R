library(ggplot2)
source('code/geom_stripes.R')
source('code/geom_effect.R')
source('code/themes.R')

TFenrichment_data=read.table('data/TF.enrichment.data.txt',header=T)
top15.mQTLs=read.table('data/top15.mQTLs.txt')
top15.eQTLs=read.table('data/top15.eQTLs.txt')
TFenrichment_data=rbind(merge(top15.mQTLs,TFenrichment_data,by.x=1,by.y=2,sort = F),merge(top15.eQTLs,TFenrichment_data,by.x=1,by.y=2,sort = F))
TFenrichment_data$top=c(rep('mQTLs\nTop enriched TFBS',30),rep('eQTLs\nTop enriched TFBS',30))
TFenrichment_data$top=factor(TFenrichment_data$top,levels=unique(TFenrichment_data$top))
colnames(TFenrichment_data)[1]='ft'
TFenrichment_data$ft=gsub('\\.1','',TFenrichment_data$ft)
TFenrichment_data$ft=factor(TFenrichment_data$ft,levels=rev(unique(TFenrichment_data$ft)))
TFenrichment_data$QTL=factor(TFenrichment_data$QTL,levels=rev(unique(TFenrichment_data$QTL)))

ggplot(data = TFenrichment_data, aes(x = estimate, y = ft)) +
    geom_effect(
        ggplot2::aes(
            xmin = ci.lb,
            xmax = ci.ub,
            color = QTL
        ),position = ggstance::position_dodgev(height = 0.5)) +
    theme_Publication()+
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("log(OR)\n(95% Confidence Interval)") +
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(0.91,0.13),legend.text = element_text(family="mono",size=12), legend.box.background = element_rect(colour = "black",size = 0.8),legend.title=element_blank())+
    scale_x_continuous(breaks = seq(-1.50,2.50,0.50),limits=c(-1.30,2.25))+
    scale_color_manual(values = c("red", "blue"))+
    facet_wrap(~top,scales='free')+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    labs(tag='d') +
    theme(plot.tag = element_text(face='bold'))
#600x450
