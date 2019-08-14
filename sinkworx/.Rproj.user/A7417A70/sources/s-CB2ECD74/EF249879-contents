# April 22, 2019

#'@title Import cell sinking data from BMG Labtech plate spectrophotometers
#'@description Transforms sinking data from plate format into tidy format. Works on .txt and .csv output from FLUOStar Omega and CLARIOStar plate readers.
#'@usage sink_bmg(data_folder, catalog = NA, time_units = "sec", RFU_cutoff = 10, skip = 1,
#'export_output = FALSE, save_as = "Sink_Database.csv")
#'@param data_folder folder containing the data files to be read, in quotations. Should be in the working directory.
#'@param catalog optional dataframe containing experimental metadata that can be merged with the imported data files. Should include a "Plate" column containing the filename (including file extension if .csv), and
#'"Well" in format A1, A2, etc.
#'@param time_units units of time in which measurements were taken (either "sec", "min", or "hr").
#'@param RFU_cutoff threshold for RFU values to be kept.
#'@param skip lines to be skipped when reading raw data files; should be sufficient to skip headers/blank lines.
#'@param export_output whether the resulting dataframe should be exported to .csv.
#'@param save_as desired filename, if results are to be exported.
#'@return A single dataframe containing tidied data from all files in data_folder, as well as experimental metadata if catalog is supplied. Creates the columns "Well", "Elapsed.Time.m", "RFU", and "Plate" (filename without extension).
#'@examples
#'sink_bmg("Sinking BMG", time_units = "min", RFU_cutoff = 7)
# --------------------------------------------------------------------------------

sink_bmg <- function(data_folder, catalog = NA, time_units = "sec", RFU_cutoff = 10, skip = 1,
                     export_output = FALSE, save_as = "Sink_Database.csv") {

# ---------------------------------------------------------------------------------

  library(magrittr)

  files <- list.files(data_folder)

  Sink_Database <- NULL

  for (f in 1:length(files)) {

    focus_file <- files[f]

    raw_data <- read.delim(file = file.path(data_folder, focus_file), sep = "\t", header = F,
                         skip = skip)

    time_intervals <- raw_data[ , 1] %>%
      .[grepl("Time", .)] %>%
      readr::parse_number(na = c(" ", "NA", "na"))

    raw_data <- raw_data %>%
      dplyr::filter(., stringr::str_detect(V1, "T") == F)

    names(raw_data) <- 1:length(raw_data)

    sink_data <- raw_data %>%
      tidyr::gather(key = "Column", value = "RFU")

    if (length(raw_data) == 3) {

      sink_data$Row <- c("A", "B")

      } else if (length(raw_data) == 6)  {

        sink_data$Row <- c("A", "B", "C", "D")

        } else if (length(raw_data) == 8) {

          sink_data$Row <- c("A", "B", "C", "D", "E", "F")

          } else if (length(raw_data) == 12) {

            sink_data$Row <- c("A", "B", "C", "D", "E", "F", "G", "H")

            } # close plate types

    sink_data <- sink_data %>%
      dplyr::mutate(., Well = paste(Row, Column, sep = ""), Row = NULL, Column = NULL) %>%
      dplyr::arrange(Well)

    if (time_units == "sec") {

      sink_data$Elapsed.Time.m <- time_intervals/60

      } else if (time_units == "min") {

        sink_data$Elapsed.Time.m <- time_intervals

        } else if (time_units == "hr") {

          sink_data$Elapsed.Time.m <- time_intervals * 60

          } # Close if statement

    sink_data$Plate <- stringr::str_remove(focus_file, ".txt")

    sink_data$RFU <- as.numeric(sink_data$RFU)

    sink_data<- sink_data %>%
      dplyr::filter(is.na(RFU) == F, RFU != "na", is.na(Elapsed.Time.m) == F, RFU >= RFU_cutoff)

    Sink_Database <- rbind(Sink_Database, sink_data)

    } # Close file loop

  if (is.na(catalog) == FALSE) {

    Sink_Database <- dplyr::inner_join(catalog, Sink_Database, by = c("Plate", "Well"))

    } # only happens if catalog is there, nothing else happens if it isn't

  if (export_output == TRUE) {

    write.csv(Sink_Database, file = save_as, row.names = FALSE)

    } # close if statement

  return(Sink_Database)

  } # end of function


#'@title Import cell sinking data from Molecular Devices Spectramax EM plate spectrophotometer
#'@description Transforms sinking data from plate format into tidy format. Works on .txt output from Spectramax EM plate reader.
#'@usage sink_spmx(data_folder, catalog = NA, time_units = "sec", RFU_cutoff = 2, skip = 3,
#'export_output = FALSE, save_as = "Sink_Database.csv" )
#'@param data_folder folder containing the data files to be read, in quotations. Should be in the working directory.
#'@param catalog optional dataframe containing experimental metadata that can be merged with the imported data files. Should include a "Plate" column containing the filename, and
#'"Well" in format A1, A2, etc.
#'@param time_units units of time in which measurements were taken (either "sec", "min", or "hr").
#'@param RFU_cutoff threshold for RFU values to be kept.
#'@param skip lines to be skipped when reading raw data files; should be sufficient to skip headers/blank lines.
#'@param export_output whether the resulting dataframe should be exported to .csv.
#'@param save_as desired filename, if results are to be exported.
#'@return A single dataframe containing tidied data from all files in data_folder, as well as experimental metadata if catalog is supplied. Creates the columns "Well", "Elapsed.Time.m", "RFU", and "Plate" (filename without extension).
#'@examples
#'sink_spmx("Spectramax data", culture_catalog, time_units = "min", RFU_cutoff = 7, save_as = "Sinking.csv")
#'
#'sinking_data <- sink_spmx("Spectramax data", export_output = FALSE)
#'
#'sinking_data <- sink_spmx("Spectramax data", culture_catalog, save_as = "Sinking_Data.csv")

# -----------------------------------------------------------------------------------

sink_spmx <- function(data_folder, catalog = NA, time_units = "sec", RFU_cutoff = 2, skip = 3,
                      export_output = FALSE, save_as = "Sink_Database.csv") {

# -----------------------------------------------------------------------------

  library(magrittr)

  files <- list.files(data_folder)

  Sink_Database <- NULL

  for (f in 1:length(files)) {

    focus_file <- files[f]

    raw_data <- read.delim(file = file.path(data_folder, focus_file), header = F, skip = skip,
                           fileEncoding = "UTF-16LE")

    time_intervals <- raw_data[ , 1] %>%
      readr::parse_time(., format = "%H:%M:%S") %>%
      na.omit(.)

    raw_data <- raw_data[ , 3:length(raw_data)]

    raw_data <- janitor::remove_empty(raw_data)

    names(raw_data) <- 1:length(raw_data)

    sink_data <- raw_data %>%
      tidyr::gather(key = "Column", value = "RFU")

    if (length(raw_data) == 3) {

      sink_data$Row <- c("A", "B")

      } else if (length(raw_data) == 6)  {

        sink_data$Row <- c("A", "B", "C", "D")

        } else if (length(raw_data) == 8) {

          sink_data$Row <- c("A", "B", "C", "D", "E", "F")

          } else if (length(raw_data) == 12) {

            sink_data$Row <- c("A", "B", "C", "D", "E", "F", "G", "H")

            } # close plate types

    sink_data <- sink_data %>%
      dplyr::mutate(., Well = paste(Row, Column, sep = ""), Row = NULL, Column = NULL) %>%
      dplyr::arrange(Well)

    sink_data$Elapsed.Time.m <- time_intervals/60

    sink_data$Plate <- stringr::str_remove(focus_file, ".txt")

    sink_data$RFU <- as.numeric(sink_data$RFU)

    sink_data<- sink_data %>%
      dplyr::filter(is.na(RFU) == F, RFU != "na", is.na(Elapsed.Time.m) == F, RFU >= RFU_cutoff)

    Sink_Database <- rbind(Sink_Database, sink_data)

  } # close file loop

  if (is.na(catalog) == FALSE) {

    Sink_Database <- dplyr::inner_join(catalog, Sink_Database, by = c("Plate", "Well"))

    } # only happens if catalog is there, nothing else happens if it isn't

  if (export_output == TRUE) {

    write.csv(Sink_Database, file = save_as, row.names = FALSE)

    } # close if statement

  return(Sink_Database)

} # end of function
