\name{article_to_df_dev}
\alias{article_to_df_dev}
\title{
Extract Valuable Data from a PubMed Record
}
\description{
Extract publication-specific information from a PubMed record driven by XML tags. The input record is a string (character-class vector of length 1) and includes PubMed-specific XML tags. Data are returned as a data frame where each row corresponds to one of the authors of the PubMed article.
}
\usage{
article_to_df_dev(pubmedArticle, autofill = FALSE, max_chars = 500, getKeywords = FALSE)
}
\arguments{
  \item{pubmedArticle}{
String including one PubMed record.
}
  \item{autofill}{
Logical. If TRUE, missing affiliations are automatically imputed based on other non-NA addresses from the same record.
}
  \item{max_chars}{
Numeric (>=0). Maximum number of characters to be extracted from the Article Abstract field.
}
  \item{getKeywords}{
Logical. If TRUE, Keywords and MeSH terms were retrieved
}
}
\details{
This function is obsolete. Please, use article_to_df() function. Right now, this is merely a wrapper for article_to_df() from easyPubMed.
}
\value{
Data frame including the extracted features. Each row correspond a different author.
}
\references{
\url{http://www.biotechworld.it/bioinf/2016/01/05/querying-pubmed-via-the-easypubmed-package-in-r/}
}
\author{
Damiano Fantini <"damiano.fantini@gmail.com">
}
\examples{
#
# Query PubMed, retrieve a selected citation and format it as a data frame
dami_query <- "Damiano Fantini[AU]"
dami_on_pubmed <- easyPubMed::get_pubmed_ids(dami_query)
dami_abstracts_xml <- easyPubMed::fetch_pubmed_data(dami_on_pubmed)
dami_abstracts_list <- easyPubMed::articles_to_list(dami_abstracts_xml)
article_to_df_dev(pubmedArticle = dami_abstracts_list[[4]], autofill = FALSE, max_chars = 100)
article_to_df_dev(pubmedArticle = dami_abstracts_list[[4]], autofill = TRUE, max_chars = 300)[1:2,]
}
