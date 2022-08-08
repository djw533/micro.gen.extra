
#' Parse interproscan gff3 file to a dataframe
#'
#' Get all domain entries for hits from interproscan into a tidy dataframe
#'
#' @param gff_file Gff file
#' @return fasta sequence as string
#' @examples
#' parse_interproscan_gff("input.iprscan.gff")

parse_interproscan_gff <- function(gff_file) {


  annotation <- read.table(text = gsub(";","\t", readLines("all_proteins.fasta.iprscan.gff")),
                           header = F,
                           sep = "\t",
                           fill = T) %>%
    #remove sequence data from bottom of file
    filter(V9 != "") %>%
    #get the length of the protein sequence
    group_by(V1) %>%
    mutate(seq_length = max(V5)) %>%
    ungroup() %>%
    #now filter to remove polypeptide entries
    filter(V3 != "polypeptide") %>%
    tidyr::pivot_longer(cols = !c(V1,V2,V3,V4,V5,V6,V7,V8,seq_length),
                        values_to = "data",
                        names_to = "delete") %>%
    #remove the old colnames
    select(!delete) %>%
    #split out into two columns:
    mutate(new_colnames = stringr::str_split_fixed(data, "=", 2)[,1],
           value = stringr::str_split_fixed(data, "=", 2)[,2]) %>%
    #remove data column
    select(!data) %>%
    #remove empty cells:
    filter(new_colnames != "" & value != "") %>%
    #pivot
    tidyr::pivot_wider(names_from = new_colnames,
                       values_from = value,
                       values_fill = NA) %>%
    rename(seqid = V1, method = V2, class = V3, start = V4, end = V5,
           score = V6, direction = V7, frame = V8)

  return(annotation)

}
