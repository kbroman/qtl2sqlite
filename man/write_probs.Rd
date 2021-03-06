% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_probs.R
\name{write_probs}
\alias{write_probs}
\title{Write genotype probabilities to database}
\usage{
write_probs(db, probs, map = NULL, quiet = TRUE)
}
\arguments{
\item{db}{Either a filename for a SQLite database, or a connection to such a database.}

\item{probs}{Genotype probabilities, as calculated by \code{\link[qtl2geno]{calc_genoprob}}.}

\item{map}{Optional map of marker positions (a list of vectors of positions).}

\item{quiet}{If FALSE, print some information about progress.}
}
\value{
None. (See Details.)
}
\description{
Write genotype probabilities to database
}
\details{
The genotype probabilities are written to a single SQLite database file. If the file exists, it is overwritten.

The data are placed in a series of tables.
\itemize{
\item \code{probs} - the probabilities, with fields \code{mat_index}, \code{marker}, \code{prob}
\item \code{markers} - marker information, with fields \code{chr},
    \code{marker}, and \code{marker_index}. If \code{map} was provided,
    there's also a field \code{pos}.
\item \code{chr} - chromosome information, with fields \code{chr}, \code{is_x_chr}, and \code{chr_index}
\item \code{ind} - the individual IDs
\item \code{geno} - genotype codes
\item \code{alleles} - the \code{"alleles"} attribute
\item \code{attributes} - the \code{"crosstype"} and \code{"alleleprobs"} attributes
}
}
\examples{
\dontrun{
library(qtl2geno)
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
map <- insert_pseudomarkers(iron, map, error_prob=0.002)
pr <- calc_genoprob(iron, map)

write_probs("iron.sqlite", pr, map)
}

}
\seealso{
\code{\link{read_probs}}
}
