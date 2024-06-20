

#' Read in annotation from gff file into a tidy dataframe
#'
#' Read in annotation from gff file into a tidy dataframe (including all fields in the decription [9th] column in gff files)
#'
#' @param gff_file Gff file
#' @return dataframe with annotation from gff file
#' @examples
#' read_gff("input.gff")



read_gff <- function(gff_file) {

  gff.data <- read.delim(file = gff_file,
                         sep = "\t",
                         header = F,
                         fill = T,
                         comment.char = "",
                         col.names = c("contig","method","feature","start","end","score","strand","frame","desc")) %>%
    filter(desc != "") %>%
    #extend the description column into multiple columns
    mutate(desc = strsplit(as.character(desc), ";")) %>%
    tidyr::unnest(desc) %>%
    #separate by the "=" character to get new col names, then pivot wider
    tidyr::separate(col = desc, sep = "=", into = c("desc","value"), extra = "merge") %>%
    tidyr::pivot_wider(names_from = "desc", values_from = "value", values_fn = list) %>%
    #mutate across all so that lists are removed and they're just characters
    mutate(across(everything(), as.character)) %>%
    #but change the start and end back to integers
    mutate(start = as.integer(start),
           end = as.integer(end))


  return(gff.data)

}
