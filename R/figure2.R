# Figure 2 uses the main models adjusted for age, sex and education.
suppressPackageStartupMessages(library(ggplot2))
plot_font<-Sys.getenv('ERP_PLOT_FONT','sans')
rs<-c(lapply(c('CAPE_TOT','SCIP_TOT','MASC_COR'),function(v)subset(cl,Outcome==v)),
 lapply(c('CNV','MMN','PINV','N1','P2'),function(m)subset(gr,ERP==m&Site=='FCz')))
figure_data<-data.frame(Outcome=c('CAPE-42 total','SCIP total','MASC correct','CNV','MMN','PINV','N1 post-S1','P2 post-S1'),
 Standardized_adjusted_difference=vapply(rs,function(r)r$beta,0),CI95_low=vapply(rs,function(r)r$low,0),CI95_high=vapply(rs,function(r)r$high,0),p_FDR=vapply(rs,function(r)r$q,0),
 y=c(8.59,7.69,6.79,4.96,4.06,3.16,1.31,.41))
figure_data$display_label<-figure_data$Outcome
figure_data$q_label<-ifelse(figure_data$p_FDR<.001,'<.001',sub('^0','',formatC(figure_data$p_FDR,format='f',digits=3)))
# Layout refers to d; keep the analysis data unchanged outside this environment.
local({d<-figure_data;
# Layout preserved from the manuscript figure.
xmin <- -2.35; xmax <- 1.95; xq <- 1.60
p <- ggplot(d)
# Original panel coordinates shifted upward; new C contains two components.
headers <- data.frame(y=c(9.50,5.87,2.23),
 label=c('A. Clinical and cognitive measures','B. Primary ERP components (FCz)',
         'C. Exploratory ERP components (FCz)'))
for(bounds in list(c(6.29,9.07),c(2.67,5.43),c(-.07,1.79))) {
 p <- p+annotate('segment',x=0,xend=0,y=bounds[1],yend=bounds[2],
                 linewidth=.42,linetype='dashed',color='grey40')
}
for(y in c(9.20,9.80,5.56,6.18,1.92,2.54)) {
 p <- p+annotate('segment',x=xmin,xend=xmax,y=y,yend=y,linewidth=.30,color='black')
}
p <- p+
 geom_errorbar(aes(x=Standardized_adjusted_difference,xmin=CI95_low,xmax=CI95_high,y=y),
               orientation='y',width=0,linewidth=.52,color='black')+
 geom_point(aes(x=Standardized_adjusted_difference,y=y),size=2.30,shape=21,
            stroke=.55,fill='black',color='black')+
 geom_text(aes(x=xq,y=y,label=q_label),hjust=.5,size=3.20,family=plot_font,color='grey20')+
 geom_text(data=headers,aes(x=0,y=y,label=label),inherit.aes=FALSE,
           hjust=.5,size=3.70,family=plot_font)
for(y in headers$y) {
 p<-p+annotate('text',x=xq,y=y,label='italic(P)[plain(FDR)]',parse=TRUE,
               hjust=.5,size=3.65,family=plot_font)
}
p <- p+
 scale_x_continuous(breaks=c(-2,-1,0,1),limits=c(xmin,xmax),expand=expansion(mult=c(0,0)))+
 scale_y_continuous(breaks=d$y,labels=d$display_label,limits=c(-.20,10.02),expand=expansion(mult=c(0,0)))+
 labs(x='Adjusted FEP–HC difference (SD units)',y=NULL)+
 coord_cartesian(clip='off')+
 theme_classic(base_size=11,base_family=plot_font)+
 theme(axis.text.x=element_text(color='black',size=10),
       axis.text.y=element_text(color='black',size=10,hjust=1,margin=margin(r=7)),
       axis.ticks.y=element_blank(),axis.line.y=element_blank(),
       axis.line.x=element_line(linewidth=.45,color='black'),
       axis.ticks.x=element_line(linewidth=.40,color='black'),
       axis.title.x=element_text(color='black',size=10.5,margin=margin(t=7)),
       plot.margin=margin(t=8,r=8,b=7,l=8))

 ggsave(file.path(out,'Figure_2.png'),p,width=7.1,height=4.85,units='in',dpi=600,bg='white')
 ggsave(file.path(out,'Figure_2.tiff'),p,width=7.1,height=4.85,units='in',dpi=300,compression='lzw',bg='white')
 if(capabilities('cairo'))ggsave(file.path(out,'Figure_2.pdf'),p,device=grDevices::cairo_pdf,width=7.1,height=4.85,units='in',bg='white')})
write.csv(figure_data,file.path(out,'Figure_2_data.csv'),row.names=FALSE)
