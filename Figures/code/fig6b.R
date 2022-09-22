library(ggplot2)
source('code/geom_stripes.R')
source('code/geom_effect.R')
source('code/themes.R')

CRE_enrichment_data=read.table('data/cCRE.enrichments.final.txt',header=T)
gene_body_enrichment_data=read.table('data/gene_annots.enrichments.final.txt',header=T)
enrichment_data=rbind(CRE_enrichment_data,gene_body_enrichment_data)
enrichment_data=subset(enrichment_data,coloc_feature%in%c('mQTL_specific_regions','mQTL_shared_regions'))
enrichment_data$coloc_feature=gsub('mQTL_shared_regions','e/mQTL-shared',gsub('mQTL_specific_regions','mQTL-specific',enrichment_data$coloc_feature))
enrichment_data=subset(enrichment_data,feature%in%c('CTCF','enhD','enhP','prom','3utr','5utr','cds','exons','introns'))
enrichment_data$ft=unlist(lapply(c('Insulator','Distal enhancer','Proximal enhancer','Promoter','3\' UTR','5\' UTR','CDS','Exon','Intron'),function(x) rep(x,2)))
enrichment_data$ft=factor(enrichment_data$ft)
enrichment_data$class=c(rep('Gene\nregulatory region',8),rep('Gene\nbody region',10))
enrichment_data$OR_lower[log(enrichment_data$OR_lower)< (-2.5)]=exp(1)^(-2.5) # for display prurposes

ggplot(data = enrichment_data, aes(x = log(OR), y = ft)) +
    geom_effect(
        ggplot2::aes(
            xmin = log(OR_lower),
            xmax = log(OR_upper),
            color = coloc_feature
        ),position = ggstance::position_dodgev(height = 0.5)) +
    theme_Publication()+
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("log(OR)\n(95% Confidence Interval)") +
    ylab(NULL)+
    theme(axis.text.x = element_text(angle = 45, size=14, hjust=1))+
    theme(legend.position=c(1.42,0.99),legend.title = element_text(size=11),legend.box.background = element_rect(colour = "black",size=0.3), legend.text = element_text(family="mono",size=11), legend.background = element_blank())+
    scale_x_continuous(breaks = seq(-2.50,2.50,1),limits=c(-2.50,2.50))+
    scale_color_manual(name="Trait-linked mCpGs",values = c("red", "#7F007F"))+
    facet_grid(class~.,scales='free')+
    theme(strip.background = element_rect(fill='white',color='white'),strip.text = element_text(size = 12, face = "plain")) +
    theme(plot.tag = element_text(face='bold'))+theme(strip.text.y = element_text(angle = 0)) +
    labs(tag='C') +
    theme(plot.tag = element_text(face='bold'))
#450x350
