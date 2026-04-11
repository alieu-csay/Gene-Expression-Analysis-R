\name{hgu133plus2ORGANISM}
\alias{hgu133plus2ORGANISM}
\alias{hgu133plus2ORGPKG}

\title{The Organism information for hgu133plus2}
\description{
  hgu133plus2ORGANISM is an R object that contains a single item: a
  character string that names the organism for which hgu133plus2 was
  built.
  hgu133plus2ORGPKG is an R object that contains a chararcter vector
  with the name of the organism package that a chip package depends on
  for its gene-centric annotation.
}
\details{
  Although the package name is suggestive of the organism for which it
  was built, hgu133plus2ORGANISM provides a simple way to programmatically
  extract the organism name.
  hgu133plus2ORGPKG provides a simple way to programmatically
  extract the name of the parent organism package.  The parent organism
  package is a strict dependency for chip packages as this is where the
  gene cetric information is ultimately extracted from.  The full
  package name will always be this string plus the extension ".db".  But
  most programatic acces will not require this extension, so its more
  convenient to leave it out.
}

\seealso{
  \itemize{
    \item \code{\link[AnnotationDbi]{AnnotationDb-class}} for use of
          the \code{select()} interface.
  }
}

\examples{
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.

## Bimap interface:
hgu133plus2ORGANISM
hgu133plus2ORGPKG
}

\keyword{datasets}

