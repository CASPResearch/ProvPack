#' Loads in compositional data for analysis
#'
#' This loads data into the "compositional" data class of the "provenace" package. compositional data
#' is that which for each sample adds up to a constant e.g. 100% This includes Bulk wt%, trace element ppm
#' and bulk mineralogy etc... The "compositional" object this generates can
#' be used in any of the default functions in "provenance" as well as those defined in this
#' package.
#' The input file which this reads is a .csv in which each measurement is a row, and the columns
#' correspond to the different variables (e.g. minerals), classification (for the typing), and sample name.
#'
#' Typing: Provenance analysis is broadly about categorising data into different "groups". This is best illustrated
#' by example, see the worked examples in the vignettes.
#' This approach allows different groups to be intuitively compared by fixing predefined groups to a particular
#' colour when the samples are loaded. This colour scheme is the same for every subsequent analysis allowing
#' for easier comparison. For example, it allows samples from one particular formation to be coloured the
#' same in UPb age spectra, MDS plots or PCA plots. Different groups of samples (e.g. Formation) are specified in a seperate file
#' which is simply a spreadsheet of the specific types (e.g. Old Red Sandstone, Greenstone) with a corresponding hexadecimal colour code.
#' This has to be manually inputted only once for the entire dataset. If this is not done, the automatically
#' created colour scheme will not be consistent, and could repeat colours, which makes it more difficult to compare
#' different groupings from method to method. A useful tool to generate hexadecimal codes from colours is: http://htmlcolorcodes.com/.
#'
#' @param datafilename Path and filename of the .csv file containing the data
#' @param typing Path and filename of the .csv file containing the colour scheme/typing. This should have a column "Colour" and "Category". The Colour
#' column should contain hexadec codes of colours, and the Category should contain the identifier of the group
#' that each colour should correspond to. These identifiers should be contained in the input file under the column header
#' given by "type_id".
#' @param tag A string to identify the distributional object for multiproxy analyses e.g. INDSCAL
#' @param sample_id String which corresponds to the column in the input file which contains the sample names.
#' @param val_range An integer list corresponding to the columns which contain the first and last piece of compositional data.
#' @param type_id String corresponding to column which contains the different classification (e.g. Formation, Age Period, protolith)
#' @param method String indicating what the default dissimilarity measure should be. See read.compositional in "provenance"
#' @param colmap String indicating colour scheme. See read.compositional
#' @return Returns an object of type "compostional".
#' @keywords compositional
#' @seealso read.compositional, maketernaryplot
#' @examples
#' \dontrun{qems <- load_compositional(datafilename = "qemscandata.csv",typing = "typing_palette.csv",tag = "Qemscan",sample_id = "SampleNo",val_range = c(8, 47),type_id = "SandType")}
#' @export
load_compositional <-  function(datafilename,
                                tag,
                                sample_id,
                                val_range,
                                type_id,
                                typing = NULL,
                                method = NULL,
                                colmap = "rainbow") {
  comp <- list()
  class(comp) <- "compositional"
  comp$name <- tag
  table <- read.csv(datafilename, header = TRUE)
  vals <- table[val_range[[1]]:val_range[[2]]]
  rownames(vals) <- make.names(table[, sample_id])
  comp$x <- vals
  if (is.null(method)) {
    if (any(comp$x == 0)) {
      method <- "bray"
    }
    else {
      method <- "aitchison"
    }
  }
  comp$method <- method
  comp$colmap <- colmap
  if (any(comp$x == 0) & method == "aitchison") {
    stop(
      paste(
        "This dataset contains zeros and is",
        "incompatible with the 'aitchison' distance"
      )
    )
  }

  #Creates a nested list 'type_samps' with a list of samples from a given typing
  #Sample names not necessarily syntactically valid for R so are 'cleaned' here too
  uniquetypes <- as.vector(unique(table[type_id]))
  types <- table[type_id]
  typeslist <- types[, (which(names(types) == type_id))]
  uniquetypeslist <-
    uniquetypes[, (which(names(uniquetypes) == type_id))]

  typings <- lapply(uniquetypeslist, function(x) {
    unique(table[which(typeslist == x), sample_id])
  })
  names(typings) <- uniquetypeslist
  typings <- lapply(typings, make.names)
  names(typings) <- lapply(names(typings), function(x)
    if (x == "") {
      x <- "Unspecified"
      return(x)
    } else {
      return(x)
    })
  comp$typing <- typings

  if(is.null(typing)){
  print(factor(names(comp$typing), levels = names(comp$typing)))
  print(factor(rainbow(length(names(comp$typing))), levels = rainbow(length(names(comp$typing)))))
  comp$TrueOrder <-
    factor(names(comp$typing), levels = names(comp$typing))
  comp$palette <-
    factor(rainbow(length(names(comp$typing))), levels = rainbow(length(names(comp$typing))))

    }else{
    typedf <- read.csv(typing)
    typedf2 <- data.frame()
    if (nrow(typedf) == 0) {
      missingtypes <- factor(names(typings))
      missingcols <- factor(rainbow(length(missingtypes)))
      print(
        paste0(
          "Following typing not specified in palette found in dataset: ",
          missingtypes
        )
      )
      print(
        "Automatically assigning them a random colour, and adding them to the palette file. Strongly recommend manually checking the appropriateness of this colour scheme and adjusting as necessary."
      )
      typedf2 <-
        data.frame(Category = missingtypes, Colour = missingcols)
    } else {
      missingtypes <-
        names(typings)[which(!(names(typings) %in% typedf$Category))]
      missingtypes <- factor(missingtypes)
      missingcols <- factor(rainbow(length(missingtypes)))
      if (length(missingtypes) > 0) {
        print(
          paste0(
            "Following typing not specified in palette found in dataset: ",
            missingtypes
          )
        )
        print(
          "Automatically assigning them a random colour, and adding them to the palette file. Strongly recommend manually checking the appropriateness of this colour scheme and adjusting as necessary."
        )
        typedf2 <-
          data.frame(Category = unlist(list(typedf$Category, missingtypes)),
                     Colour = unlist(list(typedf$Colour, missingcols)))
      }
      else {
        typedf2 <- typedf
      }
    }
    write.csv(typedf2, typing, row.names = FALSE)
    comp$TrueOrder <-
      factor(typedf2$Category[which(typedf2$Category %in% names(comp$typing))], levels = typedf2$Category)
    comp$palette <-
      factor(typedf2$Colour[which(typedf2$Category %in% names(comp$typing))], levels = typedf2$Colour)
  }
  comp$name <- tag
  return(comp)
}
