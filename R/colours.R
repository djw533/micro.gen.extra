
#' Create a vector of "n" random colours
#'
#' Create a vector of "n" random colours using grDevices, output either hexadecimal/named colours
#'
#' @param num Number of colours to generate
#' @param hex TRUE/FALSE for hexademical colours [Default = TRUE]
#' @return A vector of colours, either hexadecimal or named colours
#' @examples
#' random_n_colours(10, TRUE)
#' random_n_colours(num = 20, hex = FALSE)
#'
#' #Names can be set to vector afterwards (e.g. to create a named set of colours to pass to ggplot2, such as scale_fill_manual())
#' #For example:
#' colours_vector <- random_n_colours(3, TRUE) # first get a vector of three hex colours
#' names(colours_vector) <- c("first_name","second_name","third_name") # then apply the names
#'
#' Example output:
#'  first_name second_name  third_name
#'   "#CDC9A5"   "#6959CD"   "#D2691E"



random_n_colours <- function(num,hex = TRUE) {
  color = grDevices::colors()[grep('gr(a|e)y|white', grDevices::colors(), invert = T)]
  col_vector <- sample(color, num)

  #set the hex color strings (if specified) - otherwise return strings:
  if (isTRUE(hex)) {
    hex_colours <- list()
    for (col in col_vector) {
      hex_colours <- append(hex_colours,gplots::col2hex(col))
      hex_colours <- as.character(hex_colours)
    }
    return(hex_colours)
  } else {
    return(col_vector)
  }
}

#' Hexcode colours (from Glimma v1-hexcol.R)
#'
#' Check if string(s) are valid hex colour representation
#'
#' @param x the colour value(s) to check.
#' @param hash TRUE/FALSE - check for leading hash [Default = TRUE]
#'
#' @return Logical vector indicating if strings(s) are valid hex representations

is.hex <- function(x,hash = TRUE) {

  if (isTRUE(hash)) {
    isHex <- grepl("^#[[:xdigit:]]{6}$", x)
    isHexWithAlpha <- grepl("^#[[:xdigit:]]{8}$", x)
  } else {
    isHex <- grepl("^[[:xdigit:]]{6}$", x)
    isHexWithAlpha <- grepl("^[[:xdigit:]]{8}$", x)
  }

  return(isHex | isHexWithAlpha)
}

