library(ggplot2)
source('code/geom_stripes.R')
source('code/geom_effect.R')
source('code/themes.R')

VEPenrichment_data=read.table('data/meta.summary.vep.full.txt',header=T)
VEPenrichment_data$ft=rep(c('CGI','CGI Shore','CGI Shelf','Insulator','Distal enhancer','Proximal enhancer','Promoter','Splice site','5\' UTR','3\' UTR','Intron','NC transcript exon','CDS'),2)
VEPenrichment_data$ft=factor(VEPenrichment_data$ft)
VEPenrichment_data$class=rep(c(rep('CpG Island (CGI)\nregion',3),rep('Gene\nregulatory region',4),rep('Gene\nbody region',6)),2)
VEPenrichment_data$class=factor(VEPenrichment_data$class,levels=unique(VEPenrichment_data$class))
VEPenrichment_data$QTL=factor(VEPenrichment_data$QTL,levels=unique(VEPenrichment_data$QTL))

ggplot(data = VEPenrichment_data, aes(x = estimate, y = ft)) +
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
    theme(legend.position=c(0.70,0.90),legend.text = element_text(family="mono",size=12), legend.box.background = element_rect(colour = "black",size=0.8),legend.title=element_blank())+
    scale_color_manual(values = c("red", "blue"))+
    facet_grid(class~.,scales='free')+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='C') +
    theme(plot.tag = element_text(face='bold'))
#400x400
