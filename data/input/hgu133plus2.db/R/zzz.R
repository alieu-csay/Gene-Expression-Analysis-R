datacache <- new.env(hash=TRUE, parent=emptyenv())

hgu133plus2 <- function() showQCData("hgu133plus2", datacache)
hgu133plus2_dbconn <- function() dbconn(datacache)
hgu133plus2_dbfile <- function() dbfile(datacache)
hgu133plus2_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
hgu133plus2_dbInfo <- function() dbInfo(datacache)

hgu133plus2ORGANISM <- "Homo sapiens"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "hgu133plus2.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    txdb <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"ChipDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, txdb, envir=ns)
    namespaceExport(ns, dbNewname)
        
    ## Create the AnnObj instances
    ann_objs <- createAnnObjs.SchemaChoice("HUMANCHIP_DB", "hgu133plus2", "chip hgu133plus2", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("hgu133plus2.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(hgu133plus2_dbconn())
}

