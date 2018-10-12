require(ggplot2)
require(reshape2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
    args[2] <-getwd()
}
    
setwd(args[2])
file<-args[1]

name <-gsub(".txt","",file)
df<-read.table(file, sep = "\t")
names(df)<- c('Dis','Up','Down','Control', 'All','Up_P','Down_P','Control_P','All_P')

dfm <- melt(df[c(1:5)], id.var = "Dis")
ggplot(dfm,aes(x=Dis,y=value,group=variable,colour=variable))+
    
    geom_point(size=0.5)+
    theme_bw()+
    scale_y_continuous(name='Cumulative Number of Peaks Per Gene')+
    scale_x_continuous(name='Distance between TSS and Peak (bp)')+
    scale_color_manual(values=c("red", "blue","black","grey"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=element_text(size=rel(0.5)),
          axis.title.x=element_text(size=rel(0.5)),
          axis.ticks = element_line( size = 0.2),
          legend.text=element_text(colour = 'black',size=rel(0.5)),
          legend.title = element_blank(),  
          legend.position="top"
          ) 
    outfile<-paste(paste(name,'_cumulative',sep=''),"pdf",sep='.')
    ggsave(file=outfile,width=6,height=4)
  

dfm <- melt(df[c(1,6:8)], id.var = "Dis")
ggplot(dfm,aes(x=Dis,y=-log10(value),group=variable,colour=variable))+
    geom_point(size=0.1)+
    theme_bw()+
    scale_y_continuous(name='-Log10(P-value)')+
    scale_x_continuous(name='Distance between TSS and Peak (bp)')+
    scale_color_manual(values=c("red", "blue","black"))+
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=element_text(size=rel(0.5)),
          axis.title.x=element_text(size=rel(0.5)),
          axis.ticks = element_line( size = 0.2),
          legend.text=element_text(colour = 'black',size=rel(0.5)),
          legend.title = element_blank(),  
          legend.position="top"
          )
    outfile<-paste(paste(name,'_P-value',sep=''),"pdf",sep='.')
    ggsave(file=outfile,width=6,height=4)
