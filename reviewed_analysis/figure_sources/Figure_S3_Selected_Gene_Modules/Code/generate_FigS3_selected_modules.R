# Scientific scatter plots using R grid graphics, with physical dimensions in pt.
# All panels use the same coordinate scale, type sizes, dot scale and line widths.
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value=TRUE)
root <- if(length(args)) args[1] else dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
dir.create(file.path(root,"Rebuilt_Output"), recursive=TRUE, showWarnings=FALSE)
library(grid)
d <- read.csv(file.path(root,'Rebuilt_Output','FigS3_supporting_data.csv'),check.names=FALSE,na.strings=c('NA',''))
ampk <- subset(d,module=='Fly homologs of human AMPK-pathway genes')
stopifnot(nrow(ampk)==131)
W <- 612; rh <- 8.6; xlim <- c(-4.7,7); ticks <- c(-4,0,4)
cols <- c('Up in CS'='#B2182B','Down in CS'='#5AB4E5','Below DE cutoffs'='#B8B8B8','Not significant'='#B8B8B8')
# Removing outer text is handled by translation/cropping, never by scaling.
TOP_CROP <- 22
AXIS_LWD <- 0.42679131031036 / 0.75 # match FigS4 axes, in physical points
txt <- function(s,x,y,size=7,face='plain',just='left',color='black') grid.text(s,x=x,y=y-TOP_CROP,default.units='native',just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=face,col=color))
line <- function(x,y,color='black',dash='solid',width=AXIS_LWD) grid.lines(x=x,y=y-TOP_CROP,default.units='native',gp=gpar(col=color,lwd=width,lty=dash))
box <- function(x,y,w,h,fill=NA,stroke=TRUE) grid.polygon(x=c(x,x+w,x+w,x),y=c(y,y,y+h,y+h)-TOP_CROP,default.units='native',gp=gpar(fill=fill,col=if(stroke)'black' else NA,lwd=AXIS_LWD))
axes <- function(x,y,w,h) line(c(x,x,x+w),c(y,y+h,y+h))
radius <- function(p) 1+2.7*sqrt(pmin(50,pmax(0,-log10(pmax(p,.Machine$double.xmin))))/50)
dot <- function(x,y,r,color,open=FALSE) grid.circle(x=x,y=y-TOP_CROP,r=unit(r,'pt'),default.units='native',gp=gpar(fill=if(open)'white' else color,col=if(open)'#555555' else NA,lwd=AXIS_LWD))
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
body <- function(data,left,top,width,nrows=nrow(data)) {
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
        }
    }
    for(t in ticks){xx<-trans(t);line(rep(xx,2),c(top+h,top+h+2.5));txt(as.character(t),xx,top+h+10,7,just='centre')}
    invisible(top+h+16)
}
header <- function(title,letter,left,top,width) {
    if(nchar(letter))txt(letter,left,top+8.5,13,'bold')
    txt(title,left+width/2,top+8.5,9,'bold','centre')
}
panel <- function(module,letter,left,top,width) {
    dd<-d[d$module==module,]
    header(module,letter,left,top,width)
    body(dd,left,top+20,width)
}

# Display-only source-component groups. All original data fields are unchanged.
ampk <- read.csv(file.path(root,'Rebuilt_Output','AMPK_grouped_gene_set_and_expression.csv'),check.names=FALSE,na.strings=c('NA',''))
stopifnot(nrow(ampk)==131,anyDuplicated(ampk$flybase_id)==0)
layouts <- list()
for(col in 1:4) {
    dd<-ampk[ampk$display_column==col,]
    yy<-0;heads<-list();rows<-list()
    groups<-unique(dd$component_group)
    for(gi in seq_along(groups)) {
        lab<-groups[gi];gg<-dd[dd$component_group==lab,]
        heads[[length(heads)+1]]<-list(label=lab,y=yy)
        yy<-yy+15.5
        for(j in seq_len(nrow(gg))) {
            rows[[length(rows)+1]]<-list(row=gg[j,],y=yy+(j-.5)*rh)
        }
        yy<-yy+nrow(gg)*rh
        if(gi<length(groups))yy<-yy+6
    }
    layouts[[col]]<-list(heads=heads,rows=rows,height=yy)
}
column_height<-max(vapply(layouts,function(z)z$height,numeric(1)))
ampk_height<-22+column_height+16
ampk_grouped_panel <- function(top,letter='D') {
    header('Fly homologs of human AMPK-pathway genes',letter,14,top,584)
    y0<-top+22
    for(col in 1:4) {
        left<-14+(col-1)*148;ww<-140;axl<-left+58;axr<-left+ww-4
        trans<-function(v)axl+(v-xlim[1])/diff(xlim)*(axr-axl)
        axes(axl,y0,axr-axl,column_height)
        line(rep(trans(0),2),c(y0,y0+column_height),color='#8A8A8A',dash='dashed')
        layout<-layouts[[col]]
        for(hh in layout$heads) {
            yy<-y0+hh$y
            box(left,yy,ww,14,fill='white',stroke=FALSE)
            txt(hh$label,left+4,yy+7,9,'bold')
        }
        for(entry in layout$rows) {
            yy<-y0+entry$y;r<-entry$row
            line(c(axl-2.5,axl),rep(yy,2));txt(r$gene,axl-5,yy,7,'italic','right')
            if(is.na(r$log2FoldChange)) {
                txt('—',(axl+axr)/2,yy,7,just='centre',color='#777777')
            } else if(is.na(r$padj)) {
                dot(trans(r$log2FoldChange),yy,1.9,'#555555',TRUE)
            } else {
                dot(trans(r$log2FoldChange),yy,radius(r$padj),cols[r$status])
            }
        }
        for(t in ticks) {
            xx<-trans(t);line(rep(xx,2),c(y0+column_height,y0+column_height+2.5))
            txt(as.character(t),xx,y0+column_height+10,7,just='centre')
        }
    }
}
page_grouped <- function(path,full=TRUE) {
    top<-if(full)378 else 69
    after_ampk<-top+ampk_height
    lower_top<-after_ampk+20
    x_title_y<-if(full)lower_top+168 else after_ampk+16
    H<-ceiling(x_title_y+13)-TOP_CROP
    quartz(type='pdf',file=path,width=W/72,height=H/72,family='Arial')
    grid.newpage();pushViewport(viewport(xscale=c(0,W),yscale=c(H,0)))
    legend()
    if(full) {
        panel('FOXO','A',14,69,282)
        panel('DNA replication','B',316,69,282)
        panel('Spargel','C',14,243,282)
        ampk_grouped_panel(top)
        panel('Chromatin','E',14,lower_top,282)
        panel('Checkpoint','F',316,lower_top,282)
    } else ampk_grouped_panel(top,letter='')
    txt(expression(log[2]~'fold change (dMIC60-CS/dMIC60-WT)'),W/2,x_title_y,8,just='centre')
    popViewport();dev.off()
    message(basename(path),': ',W,' x ',H,' pt; labels/legends 7 pt, axis title 8 pt, headings 9 pt, panel letters 13 pt')
}
page_grouped(file.path(root,'Rebuilt_Output','FigS3_selected_gene_modules.pdf'),TRUE)
page_grouped(file.path(root,'Rebuilt_Output','AMPK_signaling_components_grouped_panel.pdf'),FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'Rebuilt_Output','R_sessionInfo.txt'))
