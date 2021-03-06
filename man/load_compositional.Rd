% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositional.R
\name{load_compositional}
\alias{load_compositional}
\title{Loads in compositional data for analysis}
\usage{
load_compositional(datafilename, tag, sample_id, val_range, type_id,
  typing = NULL, method = NULL, colmap = "rainbow")
}
\arguments{
\item{datafilename}{Path and filename of the .csv file containing the data}

\item{tag}{A string to identify the distributional object for multiproxy analyses e.g. INDSCAL}

\item{sample_id}{String which corresponds to the column in the input file which contains the sample names.}

\item{val_range}{An integer list corresponding to the columns which contain the first and last piece of compositional data.}

\item{type_id}{String corresponding to column which contains the different classification (e.g. Formation, Age Period, protolith)}

\item{typing}{Path and filename of the .csv file containing the colour scheme/typing. This should have a column "Colour" and "Category". The Colour
column should contain hexadec codes of colours, and the Category should contain the identifier of the group
that each colour should correspond to. These identifiers should be contained in the input file under the column header
given by "type_id".}

\item{method}{String indicating what the default dissimilarity measure should be. See read.compositional in "provenance"}

\item{colmap}{String indicating colour scheme. See read.compositional}
}
\value{
Returns an object of type "compostional".
}
\description{
This loads data into the "compositional" data class of the "provenace" package. compositional data
is that which for each sample adds up to a constant e.g. 100% This includes Bulk wt%, trace element ppm
and bulk mineralogy etc... The "compositional" object this generates can
be used in any of the default functions in "provenance" as well as those defined in this
package.
The input file which this reads is a .csv in which each measurement is a row, and the columns
correspond to the different variables (e.g. minerals), classification (for the typing), and sample name.
}
\details{
Typing: Provenance analysis is broadly about categorising data into different "groups". This is best illustrated
by example, see the worked examples in the vignettes.
This approach allows different groups to be intuitively compared by fixing predefined groups to a particular
colour when the samples are loaded. This colour scheme is the same for every subsequent analysis allowing
for easier comparison. For example, it allows samples from one particular formation to be coloured the
same in UPb age spectra, MDS plots or PCA plots. Different groups of samples (e.g. Formation) are specified in a seperate file
which is simply a spreadsheet of the specific types (e.g. Old Red Sandstone, Greenstone) with a corresponding hexadecimal colour code.
This has to be manually inputted only once for the entire dataset. If this is not done, the automatically
created colour scheme will not be consistent, and could repeat colours, which makes it more difficult to compare
different groupings from method to method. A useful tool to generate hexadecimal codes from colours is: http://htmlcolorcodes.com/.
}
\examples{
\dontrun{qems <- load_compositional(datafilename = "qemscandata.csv",typing = "typing_palette.csv",tag = "Qemscan",sample_id = "SampleNo",val_range = c(8, 47),type_id = "SandType")}
}
\seealso{
read.compositional, maketernaryplot
}
\keyword{compositional}
