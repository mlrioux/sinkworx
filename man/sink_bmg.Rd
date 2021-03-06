% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import functions.R
\name{sink_bmg}
\alias{sink_bmg}
\title{Import cell sinking data from BMG Labtech plate spectrophotometers}
\usage{
sink_bmg(data_folder, catalog = NA, time_units = "sec", RFU_cutoff = 10, skip = 1,
export_output = FALSE, save_as = "Sink_Database.csv")
}
\arguments{
\item{data_folder}{folder containing the data files to be read, in quotations. Should be in the working directory.}

\item{catalog}{optional dataframe containing experimental metadata that can be merged with the imported data files. Should include a "Plate" column containing the filename (including file extension if .csv), and
"Well" in format A1, A2, etc.}

\item{time_units}{units of time in which measurements were taken (either "sec", "min", or "hr").}

\item{RFU_cutoff}{threshold for RFU values to be kept.}

\item{skip}{lines to be skipped when reading raw data files; should be sufficient to skip headers/blank lines.}

\item{export_output}{whether the resulting dataframe should be exported to .csv.}

\item{save_as}{desired filename, if results are to be exported.}
}
\value{
A single dataframe containing tidied data from all files in data_folder, as well as experimental metadata if catalog is supplied. Creates the columns "Well", "Elapsed.Time.m", "RFU", and "Plate" (filename without extension).
}
\description{
Transforms sinking data from plate format into tidy format. Works on .txt and .csv output from FLUOStar Omega and CLARIOStar plate readers.
}
\examples{
sink_bmg("Sinking BMG", time_units = "min", RFU_cutoff = 7)
}
