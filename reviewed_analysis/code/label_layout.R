# Shared null backgrounds belong in captions; no repeated graph context/footer labels.
# The body retains its original physical size; the compact footer adds 10/14 points.
# Specific alleles/constructs, sample sizes, and test definitions remain in captions.
LABEL_HEADER_PT <- 10
LABEL_CONTEXT_SIZE <- 7
female_context <- function(text=NULL) {
  # Shared female/null context belongs in the caption. Preserve canvas size.
  if(is.null(text)) return(list(omit=TRUE))
  list(symbol='\u2640',text=text)
}
# Shared null-background context is stated once per figure caption.
# Keep the reserved strip so existing placements and plot sizes do not change.
label_context <- function(filename) {
  s <- tolower(basename(filename))
  if (grepl('experimental_strategy|working_model|rnaseq_workflow',s)) return(NULL)
  if (grepl('timeless|mtt|transfection|wb_anti_myc|heatmap|pca|mds|kegg|individual_genes|female|male|sleep|activity|locomotion|tem_|tmrm|program_ring|string_network|selected_gene_modules|publication|stress|fig3h',s)) return(list(omit=TRUE))
  NULL
}

drawDetails.mic60_label_layout <- function(x, recording=TRUE) {
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1,heights=grid::unit(c(x$body_height,x$header/72),'in'))))
  grid::pushViewport(grid::viewport(layout.pos.row=1))
  body <- if(inherits(x$body,'patchwork')) patchwork::patchworkGrob(x$body) else if(inherits(x$body,'ggplot')) ggplot2::ggplotGrob(x$body) else x$body
  grid::grid.draw(body)
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row=2))
  family <- if(names(grDevices::dev.cur())=='pdf') 'Helvetica' else 'Arial'
  y <- grid::unit(1,'npc')-grid::unit(1,'pt')
  if(is.list(x$context) && isTRUE(x$context$omit)) {
    # Leave this strip empty; shared background information is in the caption.
  } else if(is.list(x$context)) {
    # Separate font runs preserve both Unicode support and gene-symbol italics.
    symbol <- grid::textGrob(x$context$symbol,just=c('left','top'),
      gp=grid::gpar(fontfamily='Arial Unicode MS',fontsize=LABEL_CONTEXT_SIZE))
    body_label <- grid::textGrob(x$context$text,just=c('left','top'),
      gp=grid::gpar(fontfamily=family,fontsize=LABEL_CONTEXT_SIZE))
    sw <- grid::grobWidth(symbol)
    left <- grid::unit(.5,'npc')-(sw+grid::grobWidth(body_label))/2
    grid::grid.draw(grid::editGrob(symbol,x=left,y=y))
    grid::grid.draw(grid::editGrob(body_label,x=left+sw,y=y))
  } else {
    grid::grid.text(x$context,x=.5,y=y,just=c('center','top'),
      gp=grid::gpar(fontfamily=family,fontsize=LABEL_CONTEXT_SIZE))
  }
  grid::popViewport(2)
}
# Existing generators keep their save calls; this wrapper consistently reserves
# bottom context space for both PNG and PDF, including composite/grid plots.
ggsave <- function(filename,plot=ggplot2::last_plot(),width=NA,height=NA,units='in',...) {
  context <- label_context(filename)
  if(!is.null(context)) {
    stopifnot(units=='in',is.finite(width),is.finite(height))
    plot <- grid::grob(body=plot,context=context,body_height=height,header=LABEL_HEADER_PT,cl='mic60_label_layout')
    height <- height + LABEL_HEADER_PT/72
  }
  dots <- list(...)
  if(is.list(context) && grepl('\\.pdf$',filename,ignore.case=TRUE)) {
    # Quartz embeds the female glyph; base pdf's single-byte fonts cannot.
    dots$device <- function(filename,width,height,bg='white',...) {
      grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg)
    }
  }
  do.call(ggplot2::ggsave,c(list(filename=filename,plot=plot,width=width,height=height,units=units),dots))
}
