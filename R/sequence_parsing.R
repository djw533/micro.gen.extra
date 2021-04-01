fasta_from_single_gff <- function(gff_file, fasta_dir) {

  ##first get directory name and the file name from the input string:
  gff_file_split <- unlist(str_split(gff_file,"/"))

  filename = tail(gff_file_split,n=1)
  directory = paste(gff_file_split[1:length(gff_file_split) - 1 ],collapse ="/")
  if (directory == "") { # i.e. if the gff file in THIS directory
    directory = "." # then add a "." so that the glue paste below still allows the system to find the gff file
  }

  # find "##FASTA" in gff, then tail to the end of the file from this
  system(glue("tail -n +$( expr $(grep -n '##FASTA' {directory}/{filename} | cut -d ':' -f 1) + 1) {directory}/{filename} > {fasta_dir}/{filename}.fasta"))
  fasta_seq <- as.character(unlist(read.fasta(file = glue("{fasta_dir}/{filename}.fasta"),
                                              seqtype = "DNA",
                                              as.string = T)))
  #clean up:
  #file.remove("temp.fasta")

  ### perhaps add a check to make sure that it's a singlefasta file?

  return(fasta_seq)

}


fasta_from_gff_list <- function(gff_list,gff_names,fasta_dir) {

  #gff_list = a list of the gff_files to be used - with full / relative file path
  #gff_names = a list of names to be used for each file in gff_list, in the same order



  if (dir.exists(fasta_dir)) {
    stop(glue("{fasta_dir} already exists. Exiting."))
  } else {
    dir.create(fasta_dir)
  }

  #iterate over the set of gff files
  output_fasta_list <- map(gff_list, ~ fasta_from_single_gff(.x, fasta_dir))

  #now set the names:
  names(output_fasta_list) <- gff_names


  return(output_fasta_list)


}


protein_from_sequence <- function(fasta_sequence,start,end,strand) {


  ##pull fasta sequence out of the list
  #fasta_sequence <- as.character(fasta_list_input[fasta_name])

  gene_substring <- as.list(strsplit(substr(fasta_sequence,start,end), ""))[[1]]

  ##translate using the strand information
  if (strand == "forward") {
    prot_seq <- translate(seq = gene_substring,sens = "F")
  } else if (strand == "reverse") {
    prot_seq <- translate(seq = gene_substring,sens = "R")
  } else {
    stop("Strand must be either forward or reverse")
  }

  #now unlist the protein sequence:
  prot_seq <- str_flatten(prot_seq, collapse = "")

  return(prot_seq)

}



