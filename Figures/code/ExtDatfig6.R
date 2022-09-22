library(ggplot2)
source('code/geom_stripes.R')
source('code/geom_effect.R')
source('code/themes.R')
library(cowplot)

OCenrichment_data=read.table('data/summary.DNASeq.txt',header=T)
tissues=c('Breast','Colon','Kidney','Lung','Muscle','Ovary','Prostate','Testis')
OCenrichment_data$tissue=factor(rep(unlist(lapply(seq(1,8),function(x) rep(tissues[x],8))),2))
OCenrichment_data$feature=rep(rep(tissues,8),2)

OCenrichment_data$QTL=factor(c(rep('mQTLs',64),rep('eQTLs',64)),levels=c('mQTLs','eQTLs'))
OCenrichment_data=subset(OCenrichment_data,tissue==feature)


p1 <- ggplot(data = OCenrichment_data, aes(x = logOR, y = feature)) +
    geom_effect(
        ggplot2::aes(
            xmin = logORcil,
            xmax = logORciu,
            color = QTL
        ),position = ggstance::position_dodgev(height = 0.5)) +
    theme_Publication()+
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = 0, linetype=2) +
    scale_x_continuous(breaks = c(-0.5,0,0.5,1,1.5),limits=c(-1,2))+
    xlab("log(OR)\n(95% Confidence Interval)") +
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(0.80,0.095),legend.text = element_text(family="mono",size=12), legend.box.background = element_rect(colour = "black",size=0.8),legend.title=element_blank())+
    scale_color_manual(values = c("red", "blue"))+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='a') +
    theme(plot.tag = element_text(face='bold'))
#350x450


Chmmenrichment_data=read.table('data/meta.summary.Chmm.active_states.full.txt',header=T)

Chmmenrichment_data$ft=rep(c('Active TSS (TssA)','Flanking TssA','Upstream TssA','Downstream TssA','Strong Transcription','Weak Transcription','Genic Enhancer 1','Genic Enhancer 2','Active Enhancer 1','Active Enhancer 2','Weak Enhancer','ZNF Genes + Repeats'),2)
Chmmenrichment_data$ft=factor(Chmmenrichment_data$ft,levels=rev(unique(Chmmenrichment_data$ft)))
Chmmenrichment_data$QTL=factor(c(rep('mQTLs',12),rep('eQTLs',12)),levels=c('mQTLs','eQTLs'))

p2 <- ggplot(data = Chmmenrichment_data, aes(x = estimate, y = ft)) +
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
    scale_x_continuous(breaks = c(-0.5,0,0.5,1,1.5),limits=c(-1,2))+
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(0.85,0.095),legend.text = element_text(family="mono",size=12), legend.box.background = element_rect(colour = "black",size=0.8),legend.title=element_blank())+
    scale_color_manual(values = c("red", "blue"))+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='b') +
    theme(plot.tag = element_text(face='bold'))
#450x450

plot_grid(nrow = 1,p1,p2,rel_widths = c(1.5,2))
