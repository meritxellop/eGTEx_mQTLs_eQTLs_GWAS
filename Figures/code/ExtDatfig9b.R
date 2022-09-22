library(ggplot2)
library(cowplot)
source('code/geom_stripes.R')
source('code/geom_effect.R')
source('code/themes.R')

boot_mean_data=read.table('data/CpG_methylation.boot.means.txt',header=T)
discovery=read.table('data/discovery.set.counts.txt',header=T,comment.char = '!')
boot_mean_data$Tissue=factor(unlist(lapply(seq(1,9),function(x) rep(tissues[x],6))),levels=tissues)
boot_mean_data$Tissue=factor(boot_mean_data$Tissue)

boot_mean_data.m=subset(boot_mean_data,MP%in%'mCpGs')
boot_mean_data.m$class=gsub('QTL_shared','e/mQTL-shared',gsub('QTL_specific','mQTL-specific',boot_mean_data.m$class),boot_mean_data.m$class)
boot_mean_data.m$class=factor(boot_mean_data.m$class,levels=c('mQTL-specific','e/mQTL-shared','Tested'))
boot_mean_data.e=subset(boot_mean_data,MP%in%'eGenes')
boot_mean_data.e$class=gsub('QTL_shared','e/mQTL-shared',gsub('QTL_specific','eQTL-specific',boot_mean_data.e$class),boot_mean_data.e$class)
boot_mean_data.e$class=factor(boot_mean_data.e$class,levels=c('eQTL-specific','e/mQTL-shared','Tested'))

avg=mean(subset(boot_mean_data.m,class%in%'Tested','mean')[,1])
p1 <- ggplot(data = boot_mean_data.m, aes(y = Tissue, x = mean)) +
    geom_effect(
        ggplot2::aes(
            xmin = CIL,
            xmax = CIU,
            color = class
        ),position = ggstance::position_dodgev(height = 0.5)) +
    theme_Publication()+
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = avg, linetype=2) +
    xlab("Mean Methylation\n(95% Confidence Interval)") +
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(0.70,0.90),legend.text = element_text(size=12), legend.box.background = element_rect(colour = "black",size=0.8),legend.title=element_blank())+
    scale_color_manual(values = c("red","#7F007F", "orange"))+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='b') +
    theme(plot.tag = element_text(face='bold')) +
    coord_flip()+scale_y_discrete(labels=discovery$Tissue) +
    coord_flip()
#500x350

avg=mean(subset(boot_mean_data.e,class%in%'Tested','mean')[,1])
p2 <- ggplot(data = boot_mean_data.e, aes(y = Tissue, x = mean)) +
    geom_effect(
        ggplot2::aes(
            xmin = CIL,
            xmax = CIU,
            color = class
        ),position = ggstance::position_dodgev(height = 0.5)) +
    theme_Publication()+
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = avg, linetype=2) +
    xlab("Mean Gene Expression\n(95% Confidence Interval)") +
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(0.70,0.20),legend.text = element_text(size=12), legend.box.background = element_rect(colour = "black",size=0.8),legend.title=element_blank())+
    scale_color_manual(values = c("blue","#7F007F", "orange"))+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='') +
    theme(plot.tag = element_text(face='bold')) +
    coord_flip()+scale_y_discrete(labels=discovery$Tissue) +
    coord_flip()
#500x350

boot_mean_data.z=read.table('data/CpG_methylation.boot.means.z.txt',header=T)
unique(subset(boot_mean_data.z,MP%in%'mCpGs')[subset(boot_mean_data.z,MP%in%'mCpGs')$wpvalue<0.05,'Tissue'])
unique(subset(boot_mean_data.z,MP%in%'eGenes')[subset(boot_mean_data.z,MP%in%'eGenes')$wpvalue<0.05,'Tissue'])

plot_grid(nrow = 1,p1,p2)
