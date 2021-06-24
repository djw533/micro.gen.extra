
#' Get fasta sequence from single-contig gff3 file
#'
#' Get fasta sequence from a single-contig gff3 file and return it as a string
#'
#' @param gff_file Gff file
#' @param fasta_dir Directory to put the fasta file into [Default = "."]
#' @param clean_up TRUE/FALSE - clean up intermediate files. [Default = TRUE]
#' @return fasta sequence as string
#' @examples
#' fasta_from_single_gff("input.gff")

fasta_from_single_gff <- function(gff_file, fasta_dir = ".", clean_up = TRUE) {

  ##first get directory name and the file name from the input string:
  gff_file_split <- unlist(stringr::str_split(gff_file,"/"))

  filename = tail(gff_file_split,n=1)
  directory = paste(gff_file_split[1:length(gff_file_split) - 1 ],collapse ="/")
  if (directory == "") { # i.e. if the gff file in THIS directory
    directory = "." # then add a "." so that the glue paste below still allows the system to find the gff file
  }

  # find "##FASTA" in gff, then tail to the end of the file from this
  system(glue::glue("tail -n +$( expr $(grep -n '##FASTA' {directory}/{filename} | cut -d ':' -f 1) + 1) {directory}/{filename} > {fasta_dir}/{filename}.fasta"))


  #read using seqinr
  seqinr_data <- seqinr::read.fasta(file = glue::glue("{fasta_dir}/{filename}.fasta"),seqtype = "DNA", as.string = T)

  #check if singlefasta
  if (length(seqinr_data) > 1) {
    stop("GFF3 file must be of a single contig. Stopping.")
  }


  fasta_seq <- as.character(unlist(seqinr_data)) # convert to a string


  #clean up:
  if (isTRUE(clean_up)) {
    file.remove(glue::glue("{fasta_dir}/{filename}.fasta"))
  }


  return(fasta_seq)

}


#' Get fasta sequence from several single-contig gff3 file
#'
#' Get fasta sequence from several single-contig gff3 file and return it as a list of strings
#'
#' @param gff_list a list of the gff_files to be used - with full / relative file path
#' @param gff_names  a list of names to be used for each file in gff_list, in the same order
#' @param fasta_dir Directory to put the fasta files in to [Default = "."]
#' @param clean_up TRUE/FALSE - clean up intermediate files. [Default = TRUE]
#'
#'
#' @return List of fasta sequences
#' @examples
#' fasta_from_gff_list(List of gffs, names of the gff files)


fasta_from_gff_list <- function(gff_list,gff_names,fasta_dir = ".", clean_up = TRUE) {

  #gff_list = a list of the gff_files to be used - with full / relative file path
  #gff_names = a list of names to be used for each file in gff_list, in the same order



  # if (dir.exists(fasta_dir)) {
  #   stop(glue::glue("{fasta_dir} already exists. Exiting."))
  # } else {
  #   dir.create(fasta_dir)
  # }

  #iterate over the set of gff files
  output_fasta_list <- purrr::map(gff_list, ~ fasta_from_single_gff(.x, fasta_dir, clean_up))

  #now set the names:
  names(output_fasta_list) <- gff_names


  return(output_fasta_list)


}


#' Translate protein sequences from a corresponding fasta sequences, given start, end and strand
#'
#' Translate protein sequences from a corresponding fasta sequences, given start, end and strand. Ideal for
#' getting all protein sequences within a gggenesdf if the fasta sequences from the gff3 file has been extracted.
#'
#' @param fasta_sequence String of the fasta sequence
#' @param start Integer of the start position of the feature to be translated in the fasta_sequence
#' @param end Integer of the end position of the feature to be translated in the fasta_sequence
#' @param strand Orientation of the feature to be translated ("+", or "-")
#'
#' @return Protein sequence
#' @examples
#' protein_from_sequence(List of fasta_sequence, start, end, strand)
#'
#' #Particularly useful if wanting to get all the protein sequences in one/several gff3 files for comparison, e.g. running CDhit:
#'
#' #E.g. Taking a bunch of gff files in a directory:

protein_from_sequence <- function(fasta_sequence,start,end,strand) {


  ##pull fasta sequence out of the list
  #fasta_sequence <- as.character(fasta_list_input[fasta_name])

  gene_substring <- as.list(strsplit(substr(fasta_sequence,start,end), ""))[[1]]

  ##translate using the strand information
  if (strand == "forward") {
    prot_seq <- seqinr::translate(seq = gene_substring,sens = "F")
  } else if (strand == "reverse") {
    prot_seq <- seqinr::translate(seq = gene_substring,sens = "R")
  } else {
    stop("Strand must be either forward or reverse")
  }

  #now unlist the protein sequence:
  prot_seq <- stringr::str_flatten(prot_seq, collapse = "")

  return(prot_seq)

}



