#' MixTCRclean: Plot TCR motifs
#'
#' MixTCRclean is an R package to display TCR binding motifs for a set of TCRs provided by the user.
#' Typically, the input TCRs correspond to TCRs binding to a specific epitope, or isolated in a specific experiment.
#' MixTCRclean compares V usage, J usage, CDR3 length distribution and CDR3 sequence motifs to those expected from baseline TCR repertoire.
#'
#' @param input1 Either a .csv, .txt, .tsv file or a data.frame or a list in the MixTCRclean format.
#'    Columns should ideally include 'TRAV','TRAJ','cdr3_TRA','TRBV','TRBJ','cdr3_TRB'.
#'    If TCRs from multiple sources are included, the 'model' and 'species' columns are recommended.
#'
#' @param output.path Name of the output directory (default=NULL). Created if it doesn't exist.
#'
#' @param chain Which chain(s) to plot. Options: "AB" (default), "A", or "B".
#'
#' @param use.allele Whether to keep allele-level resolution. Default is FALSE (merge alleles).
#'
#' @param correct.gene.names Whether to auto-correct V/J names using an internal dictionary. Default is TRUE.
#'
#' @param use.mouse.strain Whether to treat mouse strain-specific genes separately. Default is FALSE.
#'
#' @param check.cdr3.mode Mode for filtering incompatible CDR3 sequences. 0 = none, 1 = default QC.
#'
#' @param start.lg, end.lg CDR3 positions or lengths to use for motif building or analysis. Defaults are 1 and 2.
#'
#' @param filename.output Custom filename for outputs. Default = NULL, filenames based on `model`.
#'
#' @param set.cdr3a.length,set.cdr3b.length Length to use for motif visualization. Default = NA (auto-detect).
#'
#' @param species.default Default species used if not in input1. Either "HomoSapiens" or "MusMusculus". Default is "HomoSapiens".
#'
#' @param verbose Verbosity level in the console. 0 = silent, 1 = basic info (default), 2 = detailed, 3 = full.
#'
#' @param build.clones Whether to reconstruct clones from clone.id or complex.id fields. Default is FALSE.
#'
#' @param keep.incomplete.chain Whether to retain incomplete alpha/beta chains. Default is TRUE.
#'
#' @param seq.protocol Sequencing protocol. Options: "Default" or "SEQTR". Affects gene merging.
#'
#' @param keep.colnames.origin Whether to preserve original column names. Default is FALSE.
#'
#' @return A list containing plots (if enabled), processed data, and model statistics.
#'
#' @export
#'
#'
MixTCRclean <- function(input1, output.path=NULL, chain="AB",
                      use.allele=F, correct.gene.names=T, use.mouse.strain=F, check.cdr3.mode=1, start.lg=1, end.lg=2,
                      filename.output=NULL,
                      set.cdr3a.length=NA, set.cdr3b.length=NA,
                      species.default="HomoSapiens", verbose=1, build.clones=F,
                      keep.incomplete.chain=T, seq.protocol="Default",
                      keep.colnames.origin=F){


  #######
  # Check input parameters
  #######


  if(!seq.protocol %in% c("Default", "SEQTR")){
    print("Invalid value for seq.protocol. Default value of \"Default\" will be used")
    seq.protocol <- "Default"
  }

  if (is.null(output.path)){
    print("No output.path, results will only be returned in an R list")

  } else if(!is.character(output.path) | length(output.path)>1){
    stop("Invalid value for output.path. Should be a single string of character indicating path where to save the plots.")
  }

  if(!is.logical(use.allele)){
    print("Invalid value for use.allele. Default value of FALSE will be used")
    use.allele <- F
  }

  if(!is.logical(correct.gene.names)){
    print("Invalid value for correct.gene.names. Default value of TRUE will be used")
    correct.gene.names <- T
  }
  if(correct.gene.names==F){
    print("Warning: the option to correct V/J names is turned off. You need to be sure all your V/J names follow the IMGT nomenclature")
  }

  if(!is.logical(use.mouse.strain)){
    print("Invalid value for use.mouse.strain. Default value of FALSE will be used")
    use.mouse.strain <- F
  }

  if(!check.cdr3.mode %in% c(0,1,2)){
    print("Invalid value for check.cdr3.mode. Default value of 1 will be used")
    check.cdr3.mode <- 1
  }

  if(!start.lg %in% 0:3){
    print("Invalid value for start.lg. Default value of 1 will be used")
    start.lg <- 1
  }

  if(!end.lg %in% 0:5){
    print("Invalid value for end.lg. Default value of 2 will be used")
    end.lg <- 2
  }
  if(!is.logical(build.clones)){
    print("Invalid value for build.clones. Default value of FALSE will be used")
    build.clones <- F
  }


  if(!is.na(set.cdr3a.length)){
    if(is.numeric(set.cdr3a.length)==F | set.cdr3a.length < Lmin | set.cdr3a.length > Lmax | set.cdr3a.length%%1 != 0){
      print(paste("Invalid value for set.cdr3a.length",set.cdr3a.length,". Default value will be used.", sep=""))
      set.cdr3a.length=NA
    }
  }
  if(!is.na(set.cdr3b.length)){
    if(is.numeric(set.cdr3b.length)==F | set.cdr3b.length < Lmin | set.cdr3b.length > Lmax | set.cdr3b.length%%1 != 0){
      print(paste("Invalid value for set.cdr3b.length=",set.cdr3b.length,". Default value will be used.", sep=""))
      set.cdr3b.length=NA
    }
  }

  if(! species.default %in% species.list){
    print("Wrong choice for species.default. Should be either \"HomoSapiens\" or \"MusMusculus\". Default value will be used")
    species.default <- "HomoSapiens"
  }

  if(! verbose %in% 0:3){
    print("Invalid value for verbose. Default value of 1 will be used")
    verbose <- 1
  }



  if(! chain %in% c("A","B","AB")){
    print("Invalid value for chain. Default value of \"AB\" will be used")
    chain <- "AB"
  }

  if(!is.null(filename.output)){
    if(!is.character(filename.output)){
      print("Invalid value for filename.output. Default value of the models will be used")
      filename.output <- NULL
    } else {
      if(grepl("\\\\|,|;| |\\*|/|\\?|#", filename.output)){
        filename.output <- gsub("\\\\|,|;| |\\*|/|\\?|#", "_", filename.output)
        print("Special characters, including /, \\, *, space,... are not supported in filename.output. Each of them will be changed into _")
        print(paste("New filename.output:", filename.output))
      }
    }
  }


  ###########
  #Set some specific values for different parameters
  ###########

  if( !is.null(output.path) ){
    if(!dir.exists(output.path) ){
      dir.create(output.path, recursive = TRUE);
    }
  }


  ############################
  # Load all the input data
  ############################

  if(is.character(input1)==T){
    filename.final <- basename(input1)
    filename.final <- sub("\\.csv$", "", filename.final)
    if(length(input1)==1){
      if(file.exists(input1)){
        ext <- tail(unlist(strsplit(input1,split=".", fixed=T)), n=1)
        if(ext=="csv"){
          input1 <- read.csv(input1)
        } else if(ext=="txt" | ext=="tsv"){
          input1 <- read.delim(input1, sep="\t", header=T)
        } else {
          stop("Invalid file format for input2")
        }
      } else {
        stop("Missing file for input1")
      }
    } else {
      stop("Invalid input1. If using a filename, it should be a single file. If using data, make sure it is a dataframe with the required fields")
    }
  } else if(is.data.frame(input1)==T){
    input1 <- as.data.frame(input1)
    #TODO: Should add condition here and make something else
    filename.final <- filename.output
    # Use as.data.frame because a tibble is also a data.frame but some of the
    # code has issues if input is a tibble instead of "simpler" data.frame (due
    # to column indexing by a single column that keep it as a tibble while it is
    # transforming it to a vector if a data.frame).
  } else {
    stop("Invalid value for input1. Should be a .csv or tab delimited .txt filename or a data.frame or a list generated by MixTCRviz")
  }

  # Check the compatibility between chain and the actual data.
  # If not compatible, try correcting chain. If impossible, stop the run
  chain <- verify.chain(input=input1, chain=chain)
  if(chain==""){
    stop("Incompatibilities between the data and chain parameter, and unable to infer the chain... check your input and the chain parameter")
  }


  if(chain=="A"){
    chain.list <- c("TRA");
    set.cdr3.length <- c(set.cdr3a.length)
  } else if(chain=="B"){
    chain.list <- c("TRB");
    set.cdr3.length <- c(set.cdr3b.length)
  } else if(chain=="AB"){
    chain.list <- c("TRA","TRB");
    set.cdr3.length <- c(set.cdr3a.length, set.cdr3b.length)
  }
  names(set.cdr3.length) <- chain.list


  # Check the input
  cat("\n####\nChecking input1:\n")
  check <- check_input(input=input1, chain = chain, name="input1",
                       species.default = species.default,
                       build.clones=build.clones)
  input1 <- check$data
  map.back.colnames <- check$col.map


  input1 <- clean_input(input=input1, use.allele=use.allele, correct.gene.names = correct.gene.names,
                        use.mouse.strain = use.mouse.strain, chain = chain, keep.incomplete.chain = keep.incomplete.chain,
                        species.default = species.default, check.cdr3.mode = check.cdr3.mode, start.lg=start.lg, end.lg=end.lg,
                        seq.protocol=seq.protocol, verbose=verbose)


  dir <- paste(output.path,"processed_data/", sep="")
  if(!dir.exists(dir)){
    dir.create(dir);
  }

  if(keep.colnames.origin){
    cn <- colnames(input1)
    if(length(map.back.colnames)>0){
      for(i in 1:length(cn)){
        if(!is.null(map.back.colnames[[cn[i]]])){
          cn[i] <- map.back.colnames[[cn[i]]]
        }
      }
    }
    colnames(input1) <- cn
  }

  # TODO: Why these lines?
  # quote_if_comma <- function(x) {  ifelse(grepl(",", x), paste0('"', x, '"'), x) }
  # input1 <- as.data.frame(lapply(input1, function(col) { if (is.character(col)) quote_if_comma(col) else col }))
  # write.csv(input1, file=paste(output.path,"/processed_data/",filename.final,"_2.csv", sep=""), quote=F, row.names = F, na = "")
  #


  if ("Comments_TRA" %in% names(input1) && is.list(input1$Comments_TRA)) {
    input1$Comments_TRA <- vapply(
      input1$Comments_TRA,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )
  }
  if ("Comments_TRB" %in% names(input1) && is.list(input1$Comments_TRB)) {
    input1$Comments_TRB <- vapply(
      input1$Comments_TRB,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )
  }
  write.csv(input1, file=paste(output.path,"/processed_data/",filename.final,".csv", sep=""), quote=T, row.names = F, na = "")


  return.list <- list()
  return.list$processed.data <- input1

  return(invisible(return.list))

}

