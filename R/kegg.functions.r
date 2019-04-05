  getNamedElement <- function (vector, name) 
  {
    if (name %in% names(vector)) return(vector[[name]])
    else return(as.character(NA))
  }
  
  download.keggrest <- function( pathway.id )
  {
    xml.fmt <- "http://rest.kegg.jp/get/%s/kgml"
    png.fmt <- "http://rest.kegg.jp/get/%s/image"
    
    xml.url <- sprintf(xml.fmt, pathway.id)
    png.url <- sprintf(png.fmt, pathway.id)
    
    pathway.img <- try({ readPNG(getURLContent(png.url)) }, silent=T )
    if( class(pathway.img)=="try-error" ) return(NULL)
    
    pathway.doc <- try({ xmlTreeParse(xml.url, getDTD = FALSE) }, silent=T )
    if( class(pathway.doc)=="try-error" ) return(NULL)
    
    pathway.doc.r <- xmlRoot(pathway.doc)
    isEntry <- sapply(xmlChildren(pathway.doc.r), xmlName) == "entry"

    # attrs <- xmlAttrs(pathway.doc.r)
    # kegg.pathwayinfo <- list(name = attrs[["name"]], title = attrs[["title"]] )

    entry=pathway.doc.r[isEntry][[10]]
    kegg.node.info <- lapply( pathway.doc.r[isEntry], function(entry)
    {
      attrs <- xmlAttrs(entry)
      entryID <- getNamedElement(attrs,"id")
      name <- getNamedElement(attrs,"name") # unname(unlist(strsplit(attrs["name"], " ")))
      type <- getNamedElement(attrs,"type")
      
      attrs <- xmlAttrs(xmlChildren(entry)$graphics)
      
      label <- strsplit(getNamedElement(attrs,"name"),", ")[[1]][1]
      label[is.na(label)] = ""
      label <- gsub("[.][.][.]", "", label)
      
      graphics <- list( name = getNamedElement(attrs,"name"), 
                        label = label,
                        x = as.integer(getNamedElement(attrs,"x")), 
                        y = as.integer(getNamedElement(attrs,"y")),  
                        type = getNamedElement(attrs,"type"), 
                        width = as.integer(getNamedElement(attrs,"width")),
                        height = as.integer(getNamedElement(attrs,"height")),        
                        fgcolor = getNamedElement(attrs,"fgcolor"),
                        bgcolor = getNamedElement(attrs,"bgcolor"))
      
      list(entryID = entryID, name = name, type = type, graphics = graphics)
    } )
    names(kegg.node.info) = paste( sapply(kegg.node.info,"[[", "name" ), sapply(kegg.node.info, function(x) x$graphics$x ), sapply(kegg.node.info, function(x) x$graphics$y ) )

    return( list( pathway.img=pathway.img,pathway.info=kegg.node.info ) )
  }