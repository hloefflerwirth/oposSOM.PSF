
plot.psf.pathway.keggrest <- function( kegg.pathway, signal.values=NULL, signal.values.lim, main="", highlight.genes=NULL )
{
	width = ncol(kegg.pathway$pathway.img)
	height = nrow(kegg.pathway$pathway.img)
	max.xy = max(width,height)

	par(mar = c(0, 0, 0, 0))
	plot(c(0-(max.xy-width)/2, width+(max.xy-width)/2), c(0-(max.xy-height)/2, height+(max.xy-height)/2), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=F)

	rasterImage(kegg.pathway$pathway.img, 0, 0, width, height, interpolate = F)



	title.node <- which( sapply( kegg.pathway$pathway.info, function(x) grepl("TITLE:", x$graphics$label ) ) )
	if(length(title.node)>0)
	{
		rect( kegg.pathway$pathway.info[[title.node]]$graphics$x-kegg.pathway$pathway.info[[title.node]]$graphics$width*0.6,
					height-kegg.pathway$pathway.info[[title.node]]$graphics$y+kegg.pathway$pathway.info[[title.node]]$graphics$height*0.6,
					kegg.pathway$pathway.info[[title.node]]$graphics$x+kegg.pathway$pathway.info[[title.node]]$graphics$width*0.6,
					height-kegg.pathway$pathway.info[[title.node]]$graphics$y-kegg.pathway$pathway.info[[title.node]]$graphics$height*0.6,
					col="gray90", lwd=0.5 )

		text( kegg.pathway$pathway.info[[title.node]]$graphics$x, height - kegg.pathway$pathway.info[[title.node]]$graphics$y, main,
					cex = .5, col = "black", family="sans" )
	}

	for( i in names( which( sapply(kegg.pathway$pathway.info,"[[", "type" ) == "gene" ) ) )
	{
		node.col = "gray90"
		if( paste( kegg.pathway$pathway.info[[i]]$name, kegg.pathway$pathway.info[[i]]$graphics$x, kegg.pathway$pathway.info[[i]]$graphics$y ) %in% highlight.genes )
		{
			node.col = "olivedrab1"

		} else if( i %in% names(signal.values) )
		{
			node.col = colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000)[999*(signal.values[i]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
		}

		rect( kegg.pathway$pathway.info[[i]]$graphics$x-kegg.pathway$pathway.info[[i]]$graphics$width*0.5,
					height-kegg.pathway$pathway.info[[i]]$graphics$y+kegg.pathway$pathway.info[[i]]$graphics$height*0.5,
					kegg.pathway$pathway.info[[i]]$graphics$x+kegg.pathway$pathway.info[[i]]$graphics$width*0.5,
					height-kegg.pathway$pathway.info[[i]]$graphics$y-kegg.pathway$pathway.info[[i]]$graphics$height*0.5,
					col=node.col, lwd=0.5 )

		text( kegg.pathway$pathway.info[[i]]$graphics$x, height - kegg.pathway$pathway.info[[i]]$graphics$y, kegg.pathway$pathway.info[[i]]$graphics$label,
					cex = ifelse( nchar(kegg.pathway$pathway.info[[i]]$graphics$label)<=6,.3,.24), col = "black", family="sans" )
	}

}


psf.signal.sheets <- function(psf.results)
{


  plot.psf.pathway.simple <- function( psf.object, signal.values, signal.values.lim, main="",
                                highlight.genes=NULL )
  {
    g <- igraph.from.graphNEL(psf.object$graph)

    nodedata <- psf.object$graph@nodeData@data
    edgedata <- psf.object$graph@edgeData@data
    node.type <- sapply( nodedata, function(x) x$kegg.type )
    edge.subtype <- sapply( edgedata, function(x) c( x$subtype1,x$subtype2) )

    edge.subtype <- lapply( edge.subtype, function(x) if("phosphorylation"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("dephosphorylation"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("ubiquitination"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("methylation"%in%x && sum(!is.na(x))==1) c(x,"activation") else x )
    edge.subtype <- lapply( edge.subtype, function(x) if("demethylation"%in%x && sum(!is.na(x))==1) c(x,"inhibition") else x )

    title.node <- grep( "TITLE:", V(g)$label )
    gene.nodes <- which(node.type%in%c("gene","ortholog"))
    comp.nodes <- which(node.type%in%c("compound"))
    process.nodes <- which( node.type%in%c("process","map") & igraph::degree(g)>0 )
    spare.nodes <- setdiff( seq(V(g)), c(title.node,gene.nodes,comp.nodes,process.nodes) )


    V(g)$names <- V(g)$label

    V(g)$x <- sapply( nodedata, function(x) as.numeric( x$kegg.gr.x ) )
    V(g)$y <- -sapply( nodedata, function(x) as.numeric( x$kegg.gr.y ) )

    V(g)$shape <- "rectangle"
    V(g)$shape[comp.nodes] <- "circle"

    V(g)$frame.color <- "black"
    V(g)$frame.color[spare.nodes] <- "gray40"

    V(g)$label.color <- "black"
    V(g)$label.color[spare.nodes] <- "gray40"

    V(g)$label[title.node] <- sub( "TITLE:", "", V(g)$label[title.node] )
    V(g)$label <- sapply( V(g)$label, function(x)
    {
      if( nchar(x) > 25 )
      {
        s <- strsplit(x," ")[[1]]
        return( paste( paste( s[ 1:ceiling(length(s)/2) ], collapse=" "  ),
                       paste( s[ (ceiling(length(s)/2)+1):length(s) ], collapse=" "  ), sep="\n" ) )

      }else return( x )
    })
    V(g)$label[title.node] = paste( V(g)$label[title.node], main, sep="\n")
#V(g)$label = 1:length(V(g))
# V(g)$label = psf.object$order$node.rank

    V(g)$label.cex <- 1
    V(g)$label.cex[title.node] <- 1.2

    V(g)$color <- "white"
    V(g)$color[title.node] <- "gray75"
    if(is.null(highlight.genes))  V(g)$color[c(gene.nodes,comp.nodes,process.nodes)] <- colorRampPalette(c("blue4","blue","gray90","orange","red4"))(1000)[999*(signal.values[c(gene.nodes,comp.nodes,process.nodes)]-signal.values.lim[1])/(signal.values.lim[2]-signal.values.lim[1])+1]
    if(!is.null(highlight.genes)) V(g)$color[which(highlight.genes)] <- "olivedrab1"


    V(g)$size <- 40
    V(g)$size2 <- 22
    V(g)$size[gene.nodes] <- 5 + 1.5 * sapply(V(g)$names[gene.nodes],nchar) # 10
    V(g)$size2[gene.nodes] <- 6
    V(g)$size[spare.nodes] <- 34
    V(g)$size2[spare.nodes] <- 12
    V(g)$size[process.nodes] <- 34
    V(g)$size2[process.nodes] <- 12
    V(g)$size[comp.nodes] <- 10

    E(g)$color <- "gray20"
    E(g)$label.color <- "black"
    E(g)$arrow.size <- 0
    E(g)$arrow.width <- 0
    E(g)$lty <- 3

    E(g)$label <- ""
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "phosphorylation" %in% x ) )] <- "+p"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "dephosphorylation" %in% x ) )] <- "-p"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "ubiquitination" %in% x ) )] <- "+u"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "methylation" %in% x ) )] <- "+m"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "demethylation" %in% x ) )] <- "-m"
    E(g)$label[which( sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) "dissociation" %in% x ) )] <- "#"



    # plot nodes only

    g2 <- delete_edges(g,E(g))
    plot( g2, asp=0, vertex.label.family="" )



    # plot edges

    V(g)$shape <- "none"
    V(g)$color <- NA
    V(g)$label <- NA

    # --> edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("activation","expression","reaction") %in% x) ) )  ]
    E(g2)$arrow.size <- 0.4
    E(g2)$arrow.width <- 0.8
    E(g2)$lty <- 1

    par(new=T)
    plot( g2, asp=0, edge.label.family="" )


    # ..> edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("indirect effect") %in% x) ) )  ]
    E(g2)$arrow.size <- 0.4
    E(g2)$arrow.width <- 0.8
    E(g2)$lty <- 2

    par(new=T)
    plot( g2, asp=0, edge.label.family="" )


    # --| edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("inhibition","repression") %in% x) ) )  ]
    E(g2)$arrow.size <- 0.1
    E(g2)$arrow.width <- 12
    E(g2)$lty <- 1

    par(new=T)
    plot( g2, asp=0, edge.label.family="" )


    # -- edges
    g2 <- g - E(g)[ which( !sapply( edge.subtype[apply( get.edgelist(g), 1, paste, collapse="|" )], function(x) any(c("binding/association","dissociation","compound") %in% x) ) )  ]
    E(g2)$arrow.size <- 0
    E(g2)$arrow.width <- 0
    E(g2)$lty <- 1

    par(new=T)
    plot( g2, asp=0, edge.label.family="" )

  }






  plot.psf.titlepage <- function( psf.object, signal.values )
  {
    layout(matrix(c(1,2,3,4,4,4,5,5,5),3,byrow=TRUE))

    par(mar=c(0,0,0,0))
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
    text(0.05, 0.94, psf.object$attrs$title , cex=2, adj=0)

    ### Population map of all nodes ###
    n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
    pw.genes <- unlist(sapply(psf.object$graph@nodeData@data,function(x)x$genes))
    pw.metagenes <- som.result$feature.BMU[names(gene.info$ids)[which(gene.info$ids %in% pw.genes)]]
    n.map[as.numeric(names(table(pw.metagenes)))] <- table(pw.metagenes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, preferences$dim.1stLvlSom)

    lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
    colr <- color.palette.portraits(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                            max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                            999 + 1]

    par(mar=c(3,6,3,6))
    image(matrix(spot.list.overexpression$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE )
    par(new=TRUE)
    plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
         xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="all genes",
         cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
    title(sub=paste("maximum =", max(n.map,na.rm=TRUE)),line=0)
    box()

    ### Population map of sink-nodes ###
    if(length(psf.object$sink.nodes)>0)
    {
      n.map <- matrix(0,preferences$dim.1stLvlSom,preferences$dim.1stLvlSom)
      pw.genes <- unlist(sapply(psf.object$graph@nodeData@data[psf.object$sink.nodes],function(x)x$genes))
      pw.metagenes <- som.result$feature.BMU[names(gene.info$ids)[which(gene.info$ids %in% pw.genes)]]
      if(length(pw.metagenes)>0)
      {
        n.map[as.numeric(names(table(pw.metagenes)))] <- table(pw.metagenes)
        n.map[which(n.map==0)] <- NA
        n.map <- matrix(n.map, preferences$dim.1stLvlSom)

        lim <- c(1,preferences$dim.1stLvlSom) + preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
        colr <- color.palette.portraits(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                                max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                                999 + 1]

        par(mar=c(3,6,3,6))
        image(matrix(spot.list.overexpression$overview.mask, preferences$dim.1stLvlSom), col="gray90", axes=FALSE )
        par(new=TRUE)
        plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
             xlab="",ylab="", xaxs="i", yaxs="i", col=colr, main="sink node genes",
             cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
        title(sub=paste("maximum =", max(n.map,na.rm=TRUE)),line=0)
        box()

      } else frame()


      ### Profile of mean sink signals ###
      ylim <- c(-2, 2)
      par(mar=c(2,7,4,5))

      barplot( sapply(signal.values,function(x) mean(log10(x$signal.at.sinks),na.rm=T) ),
                beside=TRUE, col=group.colors, names.arg=rep("",ncol(indata.ensID.m)),
                ylim=ylim, border=if (ncol(indata.ensID.m) < 80) "black" else NA )
      mtext(bquote("<log"[10] ~ "s>"), side=2, line=2.5, cex=1.5)

      ### Profile of max sink signals ###
      ylim <- c(-2, 2)
      par(mar=c(6,7,0,5))

      bar.coords <- barplot( sapply(signal.values,function(x) max(log10(x$signal.at.sinks),na.rm=T) ),
                            beside=TRUE, names.arg=rep("",ncol(indata.ensID.m)),
                            col=group.colors, ylim=ylim, border=if (ncol(indata.ensID.m) < 80) "black" else NA )
      mtext(bquote("log"[10] ~ "s"[max]), side=2, line=2.5, cex=1.5)

      if (ncol(indata.ensID.m)<100)
        text(bar.coords, par('usr')[3], labels=colnames(indata.ensID.m), srt=45, adj=c(1.1,1.1), xpd=TRUE)
    }
  }



  cat("Writing: PSF signal sheets\n|")
  for( i in 1:50 ) cat(" ");  cat("|\n|");  flush.console()

pw=10 #apop
pw=5 #bcell rec
pw=30 #TCA
pw=169 #epinet
  for(pw in 1:length(kegg.collection) )
  {

    kegg.collection[[pw]]$id = strsplit( kegg.collection[[pw]]$attrs$name, ":" )[[1]][2]
    kegg.pathway <- download.keggrest( kegg.collection[[pw]]$id )


    filename <- file.path(output.paths["PSF"], paste( make.names(names(kegg.collection)[pw]),".pdf",sep="") )
    pdf(filename, 29.7/2.54, 21/2.54)

    plot.psf.titlepage( psf.object=kegg.collection[[pw]], signal.values=psf.results[[pw]] )


    node.genes <- lapply(kegg.collection[[pw]]$graph@nodeData@data,function(x)x$genes)
    node.genes <- lapply(node.genes,function(x)names(gene.info$ids)[match(x,gene.info$ids)] )
    node.genes <- sapply(node.genes,function(x) any(!is.na(x)) )

    par(mfrow=c(1,1),mar=c(0,0,0,0))
    # plot.psf.pathway.simple( psf.object = kegg.collection[[pw]], highlight.genes=node.genes, main="genes with data" )


    node.genes <- sapply( kegg.collection[[pw]]$graph@nodeData@data[which(node.genes)], function(x) paste(x$kegg.name, x$kegg.gr.x, x$kegg.gr.y ) )

    if(!is.null(kegg.pathway))
    plot.psf.pathway.keggrest( kegg.pathway = kegg.pathway, highlight.genes=node.genes, main=paste(names(kegg.collection)[pw],"\ngenes with data" ) )

    sink.genes <- sapply(kegg.collection[[pw]]$graph@nodeData@data,function(x)x$kegg.id %in% kegg.collection[[pw]]$sink.nodes )

    par(mfrow=c(1,1),mar=c(0,0,0,0))
    # plot.psf.pathway.simple( psf.object = kegg.collection[[pw]], highlight.genes=sink.genes, main="sink nodes" )

    sink.genes <- sapply( kegg.collection[[pw]]$graph@nodeData@data[which(sink.genes)], function(x) paste(x$kegg.name, x$kegg.gr.x, x$kegg.gr.y ) )

    if(!is.null(kegg.pathway))
    plot.psf.pathway.keggrest( kegg.pathway = kegg.pathway, highlight.genes=sink.genes, main=paste(names(kegg.collection)[pw],"\nsink nodes") )



    y = list()
    x = list()
    for( m in 1:ncol(env$indata.ensID.m) )
    {
      g <- pathway.expression.mapping( kegg.collection[[pw]]$graph, fc.m = 10^env$indata.ensID.m[,m] )
      y[[m]] = psf.results[[pw]][[colnames(env$indata.ensID.m)[m]]]$signal.at.sinks
      x[[m]] = unlist( graph::nodeData(g,attr="expression") )
    }
    x = do.call(rbind,x)
    x = x[ , which( graph::degree(g)$inDegree > 0 | graph::degree(g)$outDegree > 0 ) ]
    x = x[ , which( apply( x, 2, var ) > 0 ) ]
    y = do.call(rbind,y)

    if( nrow(x)>=ncol(x) && !is.null(y) )
    {
      betas = as.matrix( lm( y ~ x )$coefficients )[-1,,drop=FALSE]
      rownames(betas) = colnames(x)

      betas <- apply( abs(betas), 1, max )
      betas <- betas[ do.call(c,graph::nodeData(g,attr="kegg.id")) ]
      names(betas) <- do.call(c,graph::nodeData(g,attr="kegg.id"))
      betas[ which(is.na(betas)) ] <- 0

      # plot.psf.pathway.simple( psf.object = kegg.collection[[pw]], signal.values = betas,
      #                   signal.values.lim = c(-1,1)*max(betas), main = "influencer genes" )

      if(!is.null(kegg.pathway))
      plot.psf.pathway.keggrest( kegg.pathway = kegg.pathway, signal.values = betas,
                                 signal.values.lim = c(-1,1)*max(betas),
                                 main=paste(names(kegg.collection)[pw],"\ninfluencer genes" ) )
    }

    for( m in 1:ncol(indata.ensID.m) )
    {
      # plot.psf.pathway.simple( psf.object = kegg.collection[[pw]], signal.values = log10( psf.results[[pw]][[m]]$signal.at.nodes ),
      #                   signal.values.lim = c(-1,1)*max( abs( log10( sapply( psf.results[[pw]], function(x) x$signal.at.nodes ) ) ) ),
      #                   main = colnames(indata.ensID.m)[m] )

      signal.values = log10( psf.results[[pw]][[m]]$signal.at.nodes )
      names(signal.values) = sapply( kegg.collection[[pw]]$graph@nodeData@data, function(x) paste(x$kegg.name, x$kegg.gr.x, x$kegg.gr.y ) )

      if(!is.null(kegg.pathway))
      plot.psf.pathway.keggrest( kegg.pathway = kegg.pathway, signal.values = signal.values,
                                 signal.values.lim = c(-1,1)*max( abs( log10( sapply( psf.results[[pw]], function(x) x$signal.at.nodes ) ) ) ),
                                 main=paste(names(kegg.collection)[pw],"\n",colnames(indata.ensID.m)[m]) )
    }

    dev.off()

    out.intervals <- round( seq( 1, length(kegg.collection), length.out=50+1 ) )[-1]
    cat( paste( rep("#",length( which( out.intervals == pw) ) ), collapse="" ) );	flush.console()
  }
  cat("|\n\n"); flush.console()


}

