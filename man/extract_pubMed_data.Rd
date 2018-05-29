\name{extract_pubMed_data}
\alias{extract_pubMed_data}
\title{
Extract Data from PubMed Records
}
\description{
Extract Data from PubMed Records
}
\usage{
extract_pubMed_data(pubMed_query, batch_size = 1000, 
                    getKeywords = FALSE, affi_regex_exclude = NULL)
}
\arguments{
  \item{pubMed_query}{
String (character vector of length 1): text of the query to be submitted to pubmed
}
  \item{batch_size}{
Numeric (integer): number of the PubMed records to be downloaded in a single iteration. Cannot exceed 5000. Recommended range: 200-2000
}
  \item{getKeywords}{
Logical. If TRUE, it forces keyword exraction.
}
  \item{affi_regex_exclude}{
String (character vector of length 1): Regex expression (case-insensitive) to be applied on the results (Address field). Records matching this expression will be filtered out. 
}
}
\value{
List including
\item{params }{List including the parameters used for running the Data Extraction run}
\item{data }{ Results of the Data extracted from the PubMed Records}
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/21/scraping-pubmed-data-via-easypubmed-xml-and-regex-in-r-for-a-targeting-campaign/}
}
\author{
Damiano Fantini}
\examples{
\dontrun{
pubMed_query <- "OGG1 AND 2013[PDAT]"
my_result <- extract_pubMed_data(my_PM_query, 50)
#timing
my_result$params$timing
#results
KEEP <- my_result$data$firstname == "Damiano"
my_result$data[KEEP, ]
#Num of emails retrieved
sum(!is.na(my_result$data$email))
}
}
