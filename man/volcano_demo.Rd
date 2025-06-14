% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_demo.R
\name{volcano_demo}
\alias{volcano_demo}
\title{volcano plot}
\usage{
volcano_demo(
  x,
  p.name = "P.Value",
  fc.name = "logFC",
  curve = T,
  p.value = 0.05,
  fc.value = 0.585,
  plot.name = NULL,
  gene.repel = NULL,
  width = 7,
  height = 4.5
)
}
\arguments{
\item{x}{deg result}

\item{p.name}{the colname of pvalue}

\item{fc.name}{the colname of fold change}

\item{p.value}{the threshold of pvalue}

\item{fc.value}{the threshold of fold change}

\item{plot.name}{if want to save the volcano plot file please assign it}

\item{gene.repel}{if want to point some special gene，please input a vector of gene symbol}
}
\value{
deg result and volcano plot
}
\description{
volcano plot
}
