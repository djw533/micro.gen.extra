

#' Run Hmmsearch
#'
#' Run hmmsearch and store tabular data frame, domain tabular output data frame, and an alignment of the domain hits in a list
#'
#' @param query fasta file or Biostrings AAStringSet
#' @param hmm_profile Hmm profile file (with at least one profile hmm)
#' @param input_format Query input format. One of ["fasta", "AAStringSet"],
#' @param output_dir Output directory.
#' @param cutoff Hmmer score cutoff. [Default = 20]
#' @param cpu Number of threads to use. [Default = 1]
#' @param clean Remove output directory and contents. [Default = TRUE]
#'
#' @return Whole sequence and per-domain hmmsearch results as dataframes, and alignment of hits as alignment in Biostrings format
#' @examples
#' hmmsearch_results <- run_hmmsearch(query = "seqfile.fasta", hmm_profile = "hmmfile.hmm)
#'
#' Access results as follows
#' hmmsearch_results$tbl
#' hmmsearch_results$dom_tbl
#' hmmsearch_results$aln


run_hmmsearch <- function(
  query,
  hmm_profile,
  input_format = c("fasta", "AAStringSet"),
  output_dir = "temp_hmmsearch",
  cutoff = 20,
  cpu = 1,
  clean = T) {

  check_alphabet <- function(string) {
    if (Biostrings::alphabet(string) == Biostrings::AA_ALPHABET) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  dir.create(output_dir)

  #write to fasta if input format is "AAStringSet"
  if (input_format == "AAStringSet") {
    #write to temp fasta file
    Biostrings::writeXStringSet(
      x = query,
      filepath = "temp.faa")
    input <- "temp.faa"
  } else if (input_format == "fasta") {
    input <- query
  } else  {
    unlink("temp_hmmsearch")
    stop("Input format must be either 'fasta' (default) or 'AAStringSet")
  }

  # now run hmmer
  # construct system call
  system_call <- c(
    "hmmsearch",
    "-o",
    glue::glue("{output_dir}/temp_results"),
    "--cpu",
    cpu,
    "-T",
    cutoff,
    "--tblout",
    glue::glue("{output_dir}/temp_results.tbl"),
    "-A",
    glue::glue("{output_dir}/temp_results.sto"),
    "--domtblout",
    glue::glue("{output_dir}/temp_results_domains.tbl"),
    hmm_profile,
    input
  )

  #now add the query and the hmm file
  system2(system_call)


  #read in and store the output
  hmmer_output <- list()

  #now read the hmmer output
  hmmer_output$tbl <- read.table(glue::glue("{output_dir}/temp_results.tbl"),
                                 header = F,
                                 fill = TRUE,
                                 stringsAsFactors = F)
  #set colnames
  colnames(hmmer_output$tbl) <- c(
    "target","accession","query","accession","fs.Evalue","fs.score","fs.bias","topdom.Evalue",
    "topdom.score","topdom.bias","exp","reg","clu","ov","env","dom","rep","inc","description of target"
  )


  hmmer_output$domtbl <- read.table(glue::glue("{output_dir}/temp_results_domains.tbl"),
                                    header = F,
                                    fill = TRUE,
                                    stringsAsFactors = F)
  colnames(hmmer_output$domtbl) <- c(
    "target","accession","tlen","query","accession","qlen","fs.Evalue","fs.score","fs.bias","dom.#","dom.of",
    "dom.cEvalue","dom.iEvalue","dom.score","dom.bias", "hmm.from","hmm.to","ali.from","ali.to","env.from","env.to",
    "acc","description of target"
  )


  system(glue::glue("sed -i '/^#=GR/ d' {output_dir}/temp_results.sto"))

  hmmer_output$aln <- Biostrings::readAAMultipleAlignment(
    filepath = glue::glue("{output_dir}/temp_results.sto"),
    format = "stockholm"
  )

  #remove temp dir recursively (if chosen)

  if (clean == TRUE) {
    unlink(output_dir,
           recursive = T)
  }

  return(hmmer_output)

}
