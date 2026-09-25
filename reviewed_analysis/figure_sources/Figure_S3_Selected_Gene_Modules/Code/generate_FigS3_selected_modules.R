# Scientific scatter plots using R grid graphics, with physical dimensions in pt.
# All panels use the same coordinate scale, type sizes, dot scale and line widths.
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value=TRUE)
root <- if(length(args)) args[1] else dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
dir.create(file.path(root,"Rebuilt_Output"), recursive=TRUE, showWarnings=FALSE)
library(grid)
d <- read.csv(file.path(root,'Rebuilt_Output','FigS3_supporting_data.csv'),check.names=FALSE,na.strings=c('NA',''))
ampk <- subset(d,module=='AMPK-associated genes')
stopifnot(nrow(ampk)==11)
W <- 612; rh <- 8.6; xlim <- c(-4.7,7); ticks <- c(-4,0,4)
cols <- c('Up in CS'='#B2182B','Down in CS'='#5AB4E5','Below DE cutoffs'='#B8B8B8','Not significant'='#B8B8B8')
# Removing outer text is handled by translation/cropping, never by scaling.
TOP_CROP <- 22
AXIS_LWD <- 0.42679131031036 / 0.75 # match FigS4 axes, in physical points
txt <- function(s,x,y,size=7,face='plain',just='left',color='black') grid.text(s,x=x,y=y-TOP_CROP,default.units='native',just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=color))
line <- function(x,y,color='black',dash='solid',width=AXIS_LWD) grid.lines(x=x,y=y-TOP_CROP,default.units='native',gp=gpar(col=color,lwd=width,lty=dash))
box <- function(x,y,w,h,fill=NA,stroke=TRUE) grid.polygon(x=c(x,x+w,x+w,x),y=c(y,y,y+h,y+h)-TOP_CROP,default.units='native',gp=gpar(fill=fill,col=if(stroke)'black' else NA,lwd=AXIS_LWD))
axes <- function(x,y,w,h) box(x,y,w,h)
radius <- function(p) 1+2.7*sqrt(pmin(50,pmax(0,-log10(pmax(p,.Machine$double.xmin))))/50)
dot <- function(x,y,r,color,open=FALSE) grid.circle(x=x,y=y-TOP_CROP,r=unit(r,'pt'),default.units='native',gp=gpar(fill=if(open)'white' else color,col=if(open)'#555555' else NA,lwd=AXIS_LWD))
# Annotate every estimable gene with BH-adjusted P < 0.05, independently of color.
# Keep full numeric precision in the CSVs; round only the displayed labels.
p_label <- function(p) {
    if(p < .001) return('adj. p < 0.001')
    paste0('adj. p = ', sub('0+$','',formatC(p,format='f',digits=5)))
}
p_annotation <- function(r,x,y,right) {
    if(is.na(r$log2FoldChange) || is.na(r$padj) || r$padj >= .05) return(invisible(NULL))
    label <- p_label(r$padj)
    label_width <- convertWidth(grobWidth(textGrob(label,gp=gpar(fontfamily='Arial',fontsize=7))), 'points', valueOnly=TRUE)
    px <- x + radius(r$padj) + 3
    if(px + label_width > right - 1)
        stop('Adjusted-P label does not fit to the right of ',r$gene)
    txt(label,px,y,7)
}
legend <- function() {
    labels <- c('Up in CS'='Higher in CS','Down in CS'='Higher in WT','Below DE cutoffs'='Below DE cutoffs')
    yy <- 35
    for(z in list(c(16,'Up in CS'),c(99,'Down in CS'),c(194,'Below DE cutoffs'))) {
        xx<-as.numeric(z[1]);dot(xx,yy,2.5,cols[z[2]]);txt(labels[z[2]],xx+7,yy)
    }
    dot(319,yy,2,'#555555',TRUE);txt('Adjusted P unavailable',326,yy)
    txt(expression(-log[10](italic(P)[adj])),14,53,7)
    for(i in seq_along(c(10,30,50))) {
        v<-c(10,30,50)[i];xx<-94+(i-1)*44;dot(xx,53,radius(10^-v),'#777777');txt(as.character(v),xx+8,53)
    }
    txt('— Fold change unavailable',240,53,7,color='#555555')
}
body <- function(data,left,top,width,nrows=nrow(data),row_height=rh) {
    rh <- row_height
    axl<-left+58;axr<-left+width-4;h<-rh*nrows
    trans<-function(v)axl+(v-xlim[1])/diff(xlim)*(axr-axl)
    axes(axl,top,axr-axl,h)
    line(rep(trans(0),2),c(top,top+h),color='#8A8A8A',dash='dashed')
    for(i in seq_len(nrow(data))) {
        yy<-top+(i-.5)*rh;r<-data[i,]
        line(c(axl-2.5,axl),rep(yy,2));txt(r$gene,axl-5,yy,7,'italic','right')
        if(is.na(r$log2FoldChange)) {
            txt('—',(axl+axr)/2,yy,7,just='centre',color='#777777')
        } else if(is.na(r$padj)) {
            dot(trans(r$log2FoldChange),yy,1.9,'#555555',TRUE)
        } else {
            dot(trans(r$log2FoldChange),yy,radius(r$padj),cols[r$status])
            p_annotation(r,trans(r$log2FoldChange),yy,axr)
        }
    }
    for(t in ticks){xx<-trans(t);line(rep(xx,2),c(top+h,top+h+2.5));txt(as.character(t),xx,top+h+10,7,just='centre')}
    invisible(top+h+16)
}
header <- function(title,letter,left,top,width) {
    if(nchar(letter))txt(letter,left,top+8.5,13,'bold')
    axl<-left+58; axr<-left+width-4
    box(axl,top,axr-axl,20,fill='#EFEFEF')
    txt(title,(axl+axr)/2,top+10,9,'bold','centre')
}
panel <- function(module,letter,left,top,width,row_height=rh) {
    dd<-d[d$module==module,]
    header(module,letter,left,top,width)
    body(dd,left,top+20,width,row_height=row_height)
}


# Membership follows the individual publication audit and AMPK_literature_selection.csv.
# Preserve the embedded figure's staggered columns and physical text/line sizes.
page_compact <- function(path, full=TRUE) {
    ww <- if(full) W else 310
    H <- if(full) 669 else 180
    quartz(type='pdf',file=path,width=ww/72,height=H/72,family='Arial')
    grid.newpage();pushViewport(viewport(xscale=c(0,ww),yscale=c(H,0)))
    if(full) {
        legend()
        # Equal column lengths and uniform compact panel gaps. Redistribute gene
        # rows vertically, preserving physical font, point, tick and line sizes.
        modules <- c('FOXO-associated genes','DNA replication/genome maintenance',
                     'Mitochondrial/metabolic genes','AMPK-associated genes',
                     'Chromatin regulation','DNA repair/checkpoints')
        column_top <- 69; column_bottom <- 650; gap <- 28
        for(column in 1:2) {
            ids <- seq(column,6,2)
            counts <- sapply(modules[ids],function(m)sum(d$module==m))
            row_height <- (column_bottom-column_top-3*20-2*gap)/sum(counts)
            stopifnot(row_height>=8.6)
            top <- column_top
            for(i in seq_along(ids)) {
                idx <- ids[i]
                panel(modules[idx],LETTERS[idx],if(column==1)14 else 316,top,282,row_height)
                top <- top+20+counts[i]*row_height+gap
            }
            stopifnot(abs(top-gap-column_bottom)<1e-8)
        }
        txt(expression(log[2]~'fold change (dMIC60-CS/dMIC60-WT)'),W/2,678,8,just='centre')
    } else {
        panel('AMPK-associated genes','D',14,30,282)
        txt(expression(log[2]~'fold change (dMIC60-CS/dMIC60-WT)'),ww/2,178,8,just='centre')
    }
    popViewport();dev.off()
}
page_compact(file.path(root,'Rebuilt_Output','FigS3_selected_gene_modules.pdf'),TRUE)
page_compact(file.path(root,'Rebuilt_Output','AMPK_literature_selected_panel.pdf'),FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'Rebuilt_Output','R_sessionInfo.txt'))
