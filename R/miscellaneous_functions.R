#' Compare two lists of lists for nested items
#'
#' Compare two lists comprised of lists to check whether elements of second-level lists are present entirely within other second-level lists.
#'
#' @param list1 First list of lists to be checked
#' @param list2 Second list of lists to be compared to
#'
#' @return list of TRUE/FALSE values stating whether lists in list1 are nested within any list in list2
#' @examples
#' nested.test.function(list1,list2)

nested.test.check <- function(list1,list2) {

  list.check <- unlist(purrr::map(list1, function(a) any(unlist(purrr::map(list2, function(b) ifelse(all(a %in% b) & a != b, T, F))))))

  return(list.check)
}


#' Parse snpEff output into tidy dataframe
#'
#' Parse snpEff output into tidy dataframe, using a snpEff vcf file as input
#'
#' @param snpEff_vcf First list of lists to be checked
#'
#' @return tidy dataframe
#' @examples
#' parse_snpEff_vcf(snpEff_vcf)



parse_snpEff_vcf <- function(snpEff_vcf) {


  vcf_df <- read.table(file = snpEff_vcf,
                           sep = "\t",
                           stringsAsFactors = F) %>%
    rename(contig = V1,
           locus = V2,
           id = V3,
           ref = V4,
           alt = V5,
           qual = V6,
           filter = V7,
           info = V8,
           format = V9,
           extra = V10) %>%
    #now select important stuff
    select(contig,locus,ref,alt,qual,info) %>%
    tidyr::separate(info,
                    sep = "[|]",
                    into = c("Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID",
                             "Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c",
                             "HGVS.p","cDNA.pos / cDNA.length","CDS.pos / CDS.length","AA.pos / AA.length",
                             "Distance","ERRORS / WARNINGS / INFO")) %>%
    #now separate out the allele:
    tidyr::separate(Allele,
                    sep = ";",
                    into = c("AB","AO","DP","QA","QR","RO","type","ANN")) %>%
    #select useful columns
    select(contig,locus,ref,alt,qual,type,Annotation,Annotation_Impact,Gene_Name,
           Feature_Type,Transcript_BioType,HGVS.p,`AA.pos / AA.length`) %>%
    #tidy up
    mutate(type = gsub("TYPE=","",type),
           Gene_Name = unlist(purrr::map(Gene_Name, ~ str_split(.x,"%2C")[[1]][2])))


  return(vcf_df)


}
