gggenes_df_from_gff_dir <- function(gff_dir) {

  gff_files <- list.files(path = gff_dir, full.names = T, pattern = ".gff")

  gggenes_df <- data.frame(matrix(nrow = 0 ,ncol = 8))

  colnames(gggenes_df) <- c("number","filename_prefix","start","end","gene","strand","direction","type")

  for (file in gff_files) {

    #set the filename
    filename = tail(unlist(stringr::str_split(file,"/")),1)
    filename_prefix = tail(unlist(stringr::str_split(filename,"[.]gff")),2)[1]

    #system(glue("python2 ~/scripts/gff_to_gggenes.py {file} temp_gggenes_input.csv"))

    #now read in this file;
    #temp_df_2 <- read.csv("temp_gggenes_input.csv")

    #read in the gff file and coerce into a format gggenes likes
    temp_df <- read.table(file = file,
                          sep = "\t",
                          comment.char = "",
                          fill = TRUE,
                          quote = "") %>%
      rename(contig = V1,method = V2,type = V3, start = V4,end = V5,strand = V6 ,direction = V7 ,score = V8, details  = V9) %>%
      filter(details != "")  %>% # remove lines where there isn't a 9th column
      mutate(filename_prefix = filename_prefix) %>%
      filter(type == "CDS") %>% # only take CDSs
      # now - if there is ONE  "ID=" string within the details column, create variable "gene" and set as the string immediately after "ID=", but before the next ";". Otherwise set as NA.
      mutate(gene = as.character(purrr::map(details, ~ ifelse(stringr::str_detect(.x, "ID=") & length(unlist(stringr::str_split(.x,"ID="))) == 2,  unlist(stringr::str_split(unlist(stringr::str_split(.x,"ID="))[2],";"))[1] , NA  )))) %>%
      # put in number for gggenes:
      tibble::rownames_to_column(var = "number") %>%
      mutate(number = as.numeric(number)) %>%
      #set the direction and strands:
      mutate(direction = ifelse(direction == "+", 1, -1)) %>%
      mutate(strand = ifelse(direction == 1, "forward","reverse")) %>%
      # now select columns to rbind in:
      select(number,filename_prefix,start,end,gene,strand,direction,type)


    gggenes_df <- rbind(gggenes_df,temp_df)

  }

  #file.remove("temp_gggenes_input.csv")


  return(gggenes_df)


}
