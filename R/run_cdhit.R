

#' Run CDhit on a list of protein sequences
#'
#' Run CDhit on a list of protein sequences, then read in the output into an R dataframe
#'
#' @param protein_sequences Named vector. List comprised of protein sequences, names are the unique identifiers for the sequences.
#' @return R dataframe of CDhit output
#' @examples
#' run_cdhit(list_of_protein_sequences)

run_cdhit <- function(protein_sequences) {

  # protein_sequences = named list of protein sequences
  #protein_sequences = input_protein_sequences

  #first write out the protein sequences temp file:

  temp_fasta_file <- "cdhit_temp_input_file.fasta"

  if (file.exists(temp_fasta_file)) {
    file.remove(temp_fasta_file)
  } else {
    file.create(temp_fasta_file)
  }

  #then go through the protein sequences file and append the sequences:
  for (i in 1:length(protein_sequences)) {
    cat(glue::glue(">{names(protein_sequences)[i]}"), append = T, file = temp_fasta_file, "\n")
    cat(glue::glue("{protein_sequences[i]}"), append = T, file = temp_fasta_file, "\n")
  }

  #now run cdhit: (obvs. need cdhit installed)
  system(glue::glue("cdhit -i {temp_fasta_file} -o cdhit_output -d 100"))

  ## now read the output in:
  cdhit_clusters_output <- read.csv("cdhit_output.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE) %>%
    #set names:
    rename(cdhit_cluster = V1, gene = V2) %>%
    #remove cluster col values if it's not a cluster (i.e. doesn't start with '>') :
    mutate(cdhit_cluster = ifelse(stringr::str_starts(cdhit_cluster,">") == TRUE, cdhit_cluster, NA),
           cdhit_cluster = gsub(">Cluster ","",cdhit_cluster)) %>%
    #now fill in the data downwards:
    tidyr::fill(cdhit_cluster, .direction = "down") %>%
    #now remove rows with empty values in "gene"
    filter(gene != "") %>%
    #now reformat the gene column:
    mutate(gene = stringr::str_split_fixed(stringr::str_split_fixed(gene, ">", 2)[,2],"[.][.][.]",2)[,1])


  #now return this dataframe
  return(cdhit_clusters_output)
}
