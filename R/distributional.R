readdistdf <-
  function(valdf,
           errdf = NULL,
           method = "KS",
           xlab = "age [Ma]",
           colmap = "rainbow") {
    # This is a modified version of the read.distributional function included in the provenance package.
    # Instead of taking a .csv it takes up to two data-frame objects and parses them into an object
    # of class "distributional" for analysis according to the provenance package.
    out <- list()
    if (method == "SH" & is.null(errdf)) {
      method <- "KS"
      print("Error file required for SH distance but none specified. Default method set to KS")

    }
    class(out) <- "distributional"
    out$method <- method
    out$x <- list()
    out$err <- list()
    out$colmap <- colmap
    dat <- valdf
    ns = length(dat)
    for (i in 1:ns) {
      out$x[[names(dat)[i]]] = dat[!is.na(dat[, i]), i]
    }
    if (!is.null(errdf)) {
      err <- errdf
      for (i in 1:ns) {
        out$err[[names(dat)[i]]] = dat[!is.na(err[, i]),
                                       i]
      }
    }
    d <- unlist(out$x)
    ng <- length(d)
    nb <- log(ng / ns, base = 2) + 1
    out$breaks <- seq(min(d), max(d), length.out = nb + 1)
    out$xlab <- xlab
    return(out)
  }
#' Loads in distributional data for analysis
#'
#' This loads data into the "distributional" data class of the "provenace" package. Distributional data
#' is that which for each sample it is the "distribution" of data which holds information. This
#' includes any "spectrum", e.g. UPb ages etc... The "distributional" object this generates can
#' be used in any of the default functions in "provenance" as well as those defined in this
#' package.
#' The input file which this reads is a .csv in which each measurement is a row, and the columns
#' correspond to the variable, classification, and sample name. Note that this is a different format
#' to the format required for the "provenance" package.
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
#' Note that this is the name of the sample (of which there are many individual measurements), NOT the id
#' of each individual measurement.
#' @param  val_id A string corresponding to the column which has the variable (e.g. Age, temperautre etc..)
#' @param type_id String corresponding to column which contains the different classification (e.g. Formation, Age Period, protolith)
#' @param err_id optional String corresponding to column which contains the 2S errors for variable. Required for Sircombe-Hazleton dissimilarity
#' @param SH optional Boolean indicating whether the default dissimilarity measure should be the Sircombe-Hazleton distance.
#' The default is FALSE indicating that the Kolmogorov-Smirnov distance is used.
#' @return Returns an object of type "distributional".
#' @keywords distributional
#' @examples
#' \dontrun{zircs <- load_distributional("zircondata.csv","typing_palette.csv","Zircon Ages","SampleNo", "PreferredAge", "SandType")}
#' \dontrun{zircs <- load_distributional("zircondata.csv","typing_palette.csv","Zircon Ages","SampleNo", "PreferredAge", "SandType", "PreferredAge_2s",SH = TRUE)}
#' @seealso read.distributional
#' @export
load_distributional <-
  function(datafilename,
           tag,
           sample_id,
           val_id,
           type_id,
           typing = NULL,
           err_id = NULL,
           SH = FALSE) {
    table <- read.csv(datafilename,
                      header = TRUE,
                      sep = ",")
    if (0 %in% table[, val_id]) {
      print(
        "WARNING: The data has 'zero' values. Check that there are no 'null' values that excel has automatically assigned to 0"
      )
    }
    #Aggregate all the individual age measurements in the spreadsheets into lists
    valagg = aggregate(table[, val_id], by = table[sample_id], FUN = "[")
    list_vals = valagg$x
    names(list_vals) <- valagg[, sample_id]
    #Turn the nested list of unequal list into equal lengths by padding with NA
    indx <- sapply(list_vals, length)
    new <- lapply(list_vals, `length<-`, max(indx))
    #Convert to data frame and write to csv
    valdf <- data.frame(new)

    if (!is.null(err_id)) {
      errsagg = aggregate(table[, err_id], by = table[sample_id], FUN =
                            "[")
      list_errs = errsagg$x
      names(list_errs) <- errsagg[, sample_id]

      indx2 <- sapply(list_errs, length)
      new2 <- lapply(list_errs, `length<-`, max(indx2))
      errordf <- data.frame(new2)

      #Uses above dataframes to generate a distributional data set 'zircs'
      if (SH) {
        meth <- "SH"
      } else {
        meth <- "KS"
      }
      dist <-
        readdistdf(valdf, method = meth, errordf)
    } else {
      if (SH) {
        meth <- "SH"
      } else {
        meth <- "KS"
      }
      dist <-
        readdistdf(valdf, method = meth)
    }

    #Creates a nested list with a list of samples from a given typing
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

    dist$typing <- typings
    if (is.null(typing)) {
      dist$TrueOrder <-
        factor(names(dist$typing), levels = names(dist$typing))
      dist$palette <-
        factor(rainbow(length(names(dist$typing))), levels = rainbow(length(names(dist$typing))))
    } else {
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


      dist$TrueOrder <-
        factor(typedf2$Category[which(typedf2$Category %in% names(dist$typing))], levels = typedf2$Category)
      dist$palette <-
        factor(typedf2$Colour[which(typedf2$Category %in% names(dist$typing))], levels = typedf2$Colour)
    }
    dist$name <- tag
    return(dist)
  }
