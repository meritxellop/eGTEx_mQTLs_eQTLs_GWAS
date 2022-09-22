library(ggplot2)
library(ggallin)
source('code/geom_stripes.R')
source('code/geom_effect.R')
source('code/themes.R')

RegElemeQTM_data=read.table('data/meta.summary.methylation.expression.logit.RegElements.txt',header=T,sep='\t')
RegElemeQTM_data$re=gsub('Enhancer','\nEnhancer',RegElemeQTM_data$re)
RegElemeQTM_data$ft=gsub('Sign','Sign (-)',RegElemeQTM_data$ft)
RegElemeQTM_data$ft=gsub('Expression','Gene Expression',RegElemeQTM_data$ft)
RegElemeQTM_data$ft=gsub('Methylation','CpG Methylation',RegElemeQTM_data$ft)
RegElemeQTM_data$re=factor(RegElemeQTM_data$re,levels=rev(c('Promoter','Proximal \nEnhancer','Insulator','Distal \nEnhancer')))

ggplot(data = RegElemeQTM_data, aes(x = estimate, y = re)) +
    geom_effect(
        ggplot2::aes(
            xmin = ci.lb,
            xmax = ci.ub,
            color = ft
        ),position = ggstance::position_dodgev(height = 0.5)) +
    theme_Publication()+
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("log(OR)\n(95% Confidence Interval)") +
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(0.78,0.17),legend.text = element_text(size=11), legend.box.background = element_rect(colour = "black",size=0.8),legend.title=element_blank())+
    #scale_x_continuous(breaks = seq(-1.50,2.50,0.50),limits=c(-1.30,2.25))+
    scale_x_continuous(trans = ssqrt_trans) +
    scale_color_manual(values = c("chartreuse3", "cornflowerblue","blueviolet","darkorange2"))+
    #facet_grid(class~.,scales='free')+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='e') +
    theme(plot.tag = element_text(face='bold'))
#425x400



