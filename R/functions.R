
#####
# Define some function
#####

#' @export
build_stat <- function(input, chain="AB", specidf="HomoSapiens", comp.VJL=0){

  # comp.VJL>=1 means we are computing length distributions and motifs knowing P(VJ)
  # It takdf some time, but still reasonable.

  chain.list <- paste("TR",strsplit(chain,"")[[1]], sep="")

  df <- list()
  df$specidf <- specidf
  if(comp.VJL==0){
    stat.list <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countCDR1", "countCDR2", "countCDR3.L", "countVJ", "countVJ.L")
  } else {
    stat.list <- c("L", "countL", "countV", "countJ", "countV.L", "countJ.L", "countL.VJ", "countCDR1", "countCDR2", "countCDR3.L", "countCDR3.VL", "countCDR3.JL", "countCDR3.VJL", "countVJ", "countVJ.L")
  }
  for(s in stat.list){
    df[[s]] <- list()
  }

  input[input==""] <- NA # Useful if using build_stat outside of MixTCRclean

  for(ch in chain.list){

    Vn <- paste(ch,"V", sep="")
    Jn <- paste(ch,"J", sep="")
    cdr3 <- paste("cdr3_",ch,sep="")

    df$countV[[ch]] <- table(input[,Vn])
    df$countJ[[ch]] <- table(input[,Jn])

    df$countL[[ch]] <- table(nchar(input[,cdr3]))
    df$L[[ch]] <- as.numeric(namdf(df$countL[[ch]]))
    if(length(df$countL[[ch]])>0){
      namdf(df$countL[[ch]]) <- paste("L",namdf(df$countL[[ch]]),sep="_")
    }
    df$countVJ[[ch]] <- table(input[,Vn], input[,Jn])

    df$countV.L[[ch]] <- list()
    df$countJ.L[[ch]] <- list()
    df$countVJ.L[[ch]] <- list()
    for(lg in df$L[[ch]]){
      ind <- which(nchar(input[,cdr3])==lg)
      lg.c <- paste("L",lg,sep="_")
      df$countV.L[[ch]][[lg.c]] <- table(input[ind,Vn])
      df$countJ.L[[ch]][[lg.c]] <- table(input[ind,Jn])
      df$countVJ.L[[ch]][[lg.c]] <- table(input[ind,Vn],input[ind,Jn])
    }

    if(comp.VJL>=1){

      for(V in namdf(df$countV[[ch]])){
        indv <- which(input[,Vn]==V)
        df$countCDR3.VL[[ch]][[V]] <- count_aa(input[indv,cdr3], keep.gap=0)
      }
      for(J in namdf(df$countJ[[ch]])){
        indj <- which(input[,Jn]==J)
        df$countCDR3.JL[[ch]][[J]] <- count_aa(input[indj,cdr3], keep.gap=0)
      }
      for(V in namdf(df$countV[[ch]])){
        indv <- which(input[,Vn]==V)

        for(J in namdf(df$countJ[[ch]])){
          indj <- which(input[indv,Jn]==J)
          ind <- indv[indj]
          s <- paste(V,J, sep="_")
          if(length(ind)>0){
            df$countL.VJ[[ch]][[s]] <- table(nchar(input[ind,cdr3]))
            if(length(df$countL.VJ[[ch]][[s]])>0){
              namdf(df$countL.VJ[[ch]][[s]]) <- paste("L", namdf(df$countL.VJ[[ch]][[s]]),sep="_")
            }
            df$countCDR3.VJL[[ch]][[s]] <- count_aa(input[ind,cdr3], keep.gap=0)
          } else {
            df$countL.VJ[[ch]][[s]] <- table(NA)
            df$countCDR3.VJL[[ch]][[s]] <- table(NA)
          }
        }
      }
    }
    if(length(df$countV[[ch]])>0){
      df$countCDR1[[ch]] <- count_aa(cdr123[[specidf]][[ch]][input[,Vn],"CDR1"], keep.gap=1)
      df$countCDR2[[ch]] <- count_aa(cdr123[[specidf]][[ch]][input[,Vn],"CDR2"], keep.gap=1)
    }
    df$countCDR3.L[[ch]] <- count_aa(input[,cdr3], keep.gap=0)


  }

  return(df)

}

#Compute the counts of each aa at each position.
#Would be better to have the option of treating gaps in CDR1/2 as separate amino acids, while missing data can be treated as 'unspecific'
count_aa <- function(cdr.seq, keep.gap=0){

  if(keep.gap == 0){    #All gaps, including "g" are treated as unspecific data. This can be useful for visualisation.
    tgap <- c(gap,"g")
    taa.list <- aa.list
  } else {   # "x" are treated as additonal aa, other 'gaps' are discarded. This is more correct for modelling CDR1/CDR2 loops.
    tgap <- c()
    taa.list <- c(aa.list,"g")
  }

  #First get the list of length
  l.seq <- nchar(cdr.seq)
  L <- sort(unique(l.seq))
  L <- L[L>0]  # This ensures that the cases of empty sequences are never counted.
  m.list <- list()
  for(lg in L){
    ind <- which(l.seq==lg)
    tcdr.seq <- cdr.seq[ind]
    m <- matrix(0, nrow=length(taa.list), ncol=lg)
    rownames(m) <- taa.list
    for(p in 1:lg){
      s <- substr(tcdr.seq,p,p)
      tb <- table(s)
      for(a in names(tb)){
        if(keep.gap==0){
          if(a %in% tgap){
            m[,p] <- m[,p] + tb[a]/20
          } else {
            m[a,p] <- m[a,p]+tb[a]
          }
        } else {  #Here "x" are treated as a separate amino acid and missing data (e.g., "*", "X",etc.) are not included
          if(a %in% gap==F){
            m[a,p] <- m[a,p]+tb[a]
          }
        }
      }
    }
    lc <- paste("L",lg,sep="_")
    m.list[[lc]] <- m
  }
  return(m.list)
}

#' @export
check_input <- function(input, chain="AB", name="input", species.default="HomoSapiens",
                        build.clones=F, verbose=1){

  #Check if some columns are missing, and add them with default values
  map.back.colnames <- list()

  chain.list <- paste("TR",strsplit(chain,split="")[[1]], sep="")

  if(is.data.frame(input)){

    format <- determine.format(input)

    #Do some corrections specific for VDJdb
    if(format=="VDJdb"){
      if("Species" %in% colnames(input)){
        colnames(input)[colnames(input)=="Species"] <- "species"
        map.back.colnames$species <- "Species"
      }
    }

    #Handle the format with alpha and beta chains on different rows (e.g., with clone.id)
    if(format %in% names(clone.format.col)){
      if(!build.clones){
        input <- stack_clones(input, format)
      } else {
        input <- merge_clones(input, format)
      }
    }

    col <- as.character(sapply(chain.list, function(x){c(paste(x,"V",sep=""), paste(x,"J",sep=""),paste("cdr3_",x,sep=""))}))

    #Deal with cases where column name do not follow the default of MixTCRclean
    #Do different corrections for single chain (e.g., allow V or CDR3) and paired data (allow only Va or CDR3a).
    for(i in 1:length(colnames(input))){
      cl <- colnames(input)[i]
      #Check if the col.names is in the mapping
      if(cl %in% names(mapping.colnames[[chain]])){
        mp <- mapping.colnames[[chain]][cl]
        #Make sure another entry does not already have the corrected name
        if( !mp %in% colnames(input)){
          map.back.colnames[[mp]] <- colnames(input)[i]
          colnames(input)[i] <- mp
        }
      }
    }

    #Check missing input
    for(cl in col){
      if(cl %in% colnames(input) == F){
        cn <- colnames(input)
        input <- cbind(input,"")
        colnames(input) <- c(cn, cl)
        print(paste("Missing",cl,"information in",name))
      }
    }
    #If the "species" column is not provided, we add a column with species.default
    #This is a bit suboptimal, but ok for now
    if("species" %in% colnames(input) == F){
      cn <- colnames(input)
      input <- cbind(input,species.default)
      colnames(input) <- c(cn, "species")
      print(paste("Using",species.default,"as species for all entries"))
    }
  } else {
    stop("Issues with input format")
  }
  lst <- list()
  lst[["data"]] <- input
  lst[["col.map"]] <- map.back.colnames
  return(lst)
}

#' @export
clean_input <- function(input, use.allele=F, correct.gene.names=T, use.mouse.strain=F,
                        chain="AB", species.default="HomoSapiens", check.cdr3.mode=2,
                        keep.incomplete.chain=T, start.lg=1, end.lg=2, seq.protocol="Default",
                        merge.ambiguous=T, verbose=1){

  ####
  # Clean the input by removing CDR3 with weird characters, longer than Lmax or shorter than Lmin
  # Correct VJ genes based on our dictionary
  # species.default is only used if input does not contain the "species" column
  # merge.ambiguous should be set to FALSE ONLY if clean_input is used outside of MixTCRclean
  # this will prevent mapping TRBV6-2 to TRBV6-2/6-3
  ####

  if(is.data.frame(input)){
    input <- as.data.frame(input)
    # Some code used here don't handle tibbles that are a type of data.frame
    # (as.data.frame drops these tibble class).
  }

  if("species" %in% colnames(input)){
    species.list <- unique(input[,"species"])
    use.species.default <- F
  } else {
    species.list <- c(species.default)
    use.species.default <- T
  }

  if(chain=="AB"){
    col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")
    segment.list <- c("TRAV","TRAJ","TRBV","TRBJ")
    cdr3.list <- c("cdr3_TRA","cdr3_TRB")
  }
  if(chain=="A"){
    col <- c("TRAV","TRAJ","cdr3_TRA")
    segment.list <- c("TRAV","TRAJ")
    cdr3.list <- c("cdr3_TRA")
  }
  if(chain=="B"){
    col <- c("TRBV","TRBJ","cdr3_TRB")
    segment.list <- c("TRBV","TRBJ")
    cdr3.list <- c("cdr3_TRB")
  }


  #Replace empty values by NA
  for(i in col){
    input[which(input[,i] == '' | input[,i] == ""),i] <- NA
  }

  #Set to NA CDR3 sequences with incompatible lengths or weird characters

  #print("Checking non-aa characters")
  for(cdr3 in cdr3.list){
    nc <- nchar(input[,cdr3])
    ind <- which( nc < Lmin | nc > Lmax | grepl('[^ACDEFGHIKLMNPQRSTVWY]', input[,cdr3]) == T)
    input[ind,cdr3] <- NA
  }

  #Remove anything that comes after parenthesis (this is the case for MiXCR data for instance)
  for(s in segment.list){
    ind <- grep("\\(",input[,s])
    input[ind,s] <- sapply(ind,function(i){strsplit(input[i,s], split="\\(")[[1]][1]})
  }

  #If multiple V or J genes (separated by "," or "\" or ";" or " or "), keep only the first one
  for(s in segment.list){
    ind <- grep("\\\\|,|;| or ",input[,s])
    if(length(ind)>0){
      cor <- sapply(ind,function(i){strsplit(input[i,s], split="\\\\|,|;| or ")[[1]][1]})
      if(verbose>=1) {
        entry <- ifelse(length(ind)==1,"entry","entries")
        print(paste(length(ind)," ",entry," with multiple ",s," segments (only the first segment will be kept):", sep=""))

        if(verbose==1 | verbose==2){
          print("Use verbose=3 to see them all")
        } else {
          mt <- cbind(input[ind,s], cor)
          colnames(mt) <- c("Original", "Corrected")
          print(mt)
        }
      }
      input[ind,s] <- cor
    }
  }
  ###############################################################################
  # Add something for the cleaning of SEQTR data when merge_TRAV is True
  # We have segments with double names thht are separated by "/" not ", " or something else
  # But we still want to keep only one name (first for e.g.) otherwise it is not corrected afterwards (by clean.name.allele, nor merge_TRAV)

  if (seq.protocol=="SEQTR" & !use.mouse.strain & species.default=="MusMusculus"){
    for(s in segment.list){
      ind <- grep("/(?!DV)", input[, s], perl = TRUE)
      if(length(ind)>0){
        cor <- sapply(ind,function(i){strsplit(input[i,s], split="/")[[1]][1]})
        if(verbose>=1) {
          entry <- ifelse(length(ind)==1,"entry","entries")
          print(paste(length(ind)," ",entry," with multiple ",s," segments separated with / (only the first segment will be kept):", sep=""))

          if(verbose==1 | verbose==2){
            print("Use verbose=3 to see them all")
          } else {
            mt <- cbind(input[ind,s], cor)
            colnames(mt) <- c("Original", "Corrected")
            print(mt)
          }
        }
        input[ind,s] <- cor
      }
    }

  }


  ###############################################################################
  #Remove spaces
  for(s in segment.list){
    input[,s] <- gsub(" ","",input[,s])
  }

  #Remove or add alleles
  if(!use.allele){
    #print("Removing alleles")
    for(s in segment.list){
      ind <- grep("*",input[,s], fixed=T)
      input[ind,s] <- gsub("\\*0[0-9]/0[0-9]", "", input[ind,s]) #This are entries with ambiguous allele assignment
      input[ind,s] <- gsub("\\*0[0-9]", "", input[ind,s])
    }
  } else{
    #print("Adding alleles")
    for(s in segment.list){
      ind <- which(!grepl("*",input[,s], fixed=T) & !is.na(input[,s]) )
      if(use.species.default){
        al <- allele.default[[seq.protocol]][[species.default]][input[ind,s]]
      } else {
        al <- sapply(ind, function(i){allele.default[[seq.protocol]][[input[i,"species"]]][input[i,s]]})
      }
      al[is.na(al)] <- "01" #This happens in case of wrong gene names, since gene names were not yet corrected
      input[ind,s] <- paste(input[ind,s], al, sep="*")
    }
  }

  if(merge.ambiguous){
    #Do a manual correction for TRBV6-2 and TRBV6-3 -> TRBV6-2/6-3
    if(chain=="B" | chain=="AB"){
      if(seq.protocol=="Default" | seq.protocol=="SEQTR"){
        ind <- which((grepl("TRBV6-2",input[,"TRBV"], fixed=T) | grepl("TRBV6-02",input[,"TRBV"], fixed=T) | grepl("TRBV6-3",input[,"TRBV"], fixed=T) | grepl("TRBV6-03",input[,"TRBV"], fixed=T))
                     & input[,"TRBV"] != "TRBV6-2/6-3" & input[,"TRBV"] != "TRBV6-2/6-3*01")
        if(length(ind)>0){
          if(verbose>0){
            print("Mapping all TRBV6-2 and TRBV6-3 to TRBV6-2/6-3, since they cannot be distinguished at the sequencing level")
          }
          if(use.allele){
            input[ind,"TRBV"] <- "TRBV6-2/6-3*01"
          } else {
            input[ind,"TRBV"] <- "TRBV6-2/6-3"
          }
        }
      }
    }

    #Do a manual correction for TRBV12-3 or TRBV12-4 -> TRBV12-3/12-4

    if(chain=="B" | chain=="AB"){
      if(seq.protocol=="SEQTR"){
        ind <- which((grepl("TRBV12-3",input[,"TRBV"], fixed=T) | grepl("TRBV12-03",input[,"TRBV"], fixed=T) | grepl("TRBV12-4",input[,"TRBV"], fixed=T) | grepl("TRBV12-04",input[,"TRBV"], fixed=T)) &
                       input[,"TRBV"] != "TRBV12-3/12-4" & input[,"TRBV"] != "TRBV12-3/12-4*01")

        if(length(ind)>0){
          if(verbose>0){
            print("Mapping all TRBV12-3 or TRBV12-4 to TRBV12-3/12-4, since they cannot be distinguished with SEQTR protocol")
          }
          if(use.allele){
            input[ind,"TRBV"] <- "TRBV12-3/12-4*01"
          } else {
            input[ind,"TRBV"] <- "TRBV12-3/12-4"
          }
        }
      } else {
        ind <- grep("TRBV12-3/12-4",input[,"TRBV"])
        if(length(ind)>0){
          if(verbose>0){
            print("*** WARNING: TCRs contain TRBV12-3/12-4 entries.")
            print("    If you are using data generated with SEQTR protocol, make sure to specify it with seq.protocol=\"SEQTR\"")
          }
        }
      }
    }
  }

  ###################
  # Correct gene names
  # If alleles, it will correct the gene name, and keep the allele. If the allele cannot be found, it will remove it
  # If genes, it will correct the gene name
  # If gene name cannot be corrected, it gives NA
  ###################

  if(correct.gene.names){
    #print("Check V/J names")
    input <- correct.VJnames(input=input, segment.list=segment.list, species.default=species.default,
                             use.allele=use.allele, seq.protocol=seq.protocol,verbose=verbose)
  } else {
    for(species in species.list){

      if(!use.species.default){
        ind.species <- which(input[,"species"]==species)
      } else {
        ind.species <- 1:dim(input)[1]
      }
      for(s in segment.list){
        if(use.allele){
          name.list <- gene.allele.list[[seq.protocol]][[species]][substr(gene.allele.list[[seq.protocol]][[species]],1,4)==s]
        } else {
          name.list <- gene.list[[seq.protocol]][[species]][substr(gene.list[[seq.protocol]][[species]],1,4)==s]
        }
        ind <- which(input[ind.species,s] %in% name.list == F & !is.na(input[ind.species,s]) )
        if(length(ind)>=1 & verbose != 0){
          missing <- sort(unique(input[ind.species[ind],s]))
          nm <- ifelse(length(missing)==1,"name","names")
          cat("\n")
          print(paste("*** ",length(missing)," ",s," ",nm," in ",length(ind)," entries absent from IMGT in ",species," ***",sep=""))
          print(missing)
          cat("\n")
        }
        input[ind.species[ind],s] <- NA
      }
    }
  }

  if(check.cdr3.mode > 0){
    #print("Checking CDR3")
    input <- check_cdr3(input=input, chain=chain, species.default=species.default, check.cdr3.mode=check.cdr3.mode, start.lg=start.lg, end.lg=end.lg, verbose=verbose)
  }

  ################
  # Do an extra correction for mouse entries, where only gene level analyses are allowed
  # and TRAV genes can be merged
  ################

  if(!use.species.default){
    ind <- which(input[,"species"]=="MusMusculus")
  } else {
    if(species.default=="MusMusculus"){
      ind <- 1:dim(input)[1]
    } else {
      ind <- c()
    }
  }
  if(length(ind)>0){

    if(use.allele){
      #Remove the alleles (if(use.allele==F), this was done before)
      for(s in segment.list){
        input[ind,s] <- unlist(lapply(input[ind,s], function(x){unlist(strsplit(x,split="*", fixed=T))[1]}))
      }
    }

    if(!use.mouse.strain & chain != "B"){
      input[ind,] <- merge_mouse_TRAV(input[ind,])  #WARNING: This only works if alleles have been removed (so far always the case in mouse)
      }
  }

  if(!keep.incomplete.chain){
    chain.list <- paste("TR",strsplit(chain,split="")[[1]], sep="")
    for(ch in chain.list){
      cl <- c(paste(ch,"V",sep=""),paste(ch,"J",sep=""),paste("cdr3_",ch,sep=""))
      ind <- apply(input[,cl],1,function(x){any(is.na(x))})
      input[ind,cl] <- NA
    }
  }

  return(input)

}

#' @export
check_cdr3 <- function(input, chain="AB", species.default="HomoSapiens", check.cdr3.mode=2, start.lg=1, end.lg=2, verbose=1){

  # Clean the CDR3 based on the V and J usage.
  # This should be applied after correcting the gene names, and adding the species if needed
  # species.default is only used if es.all does not contain the "species" column
  # If the allele is given in the gene name, the allele will be used.

  use.species.default <- F
  if("species" %in% colnames(input)){
    species.list <- unique(input[,"species"])
  } else {
    species.list <- c(species.default)
    use.species.default <- T
  }

  chain.list <- paste("TR",strsplit(chain,split="")[[1]], sep="")


  for(ch in chain.list){

    V <- paste(ch,"V",sep="")
    J <- paste(ch,"J",sep="")
    cdr3 <- paste("cdr3_",ch,sep="")

    for(species in species.list){

      if(!use.species.default){
        ind.species <- which(input[,"species"]==species)
      } else{
        ind.species <- 1:dim(input)[1]
      }

      if(check.cdr3.mode==0){
        ind.first <- c()
        ind.last <- c()
      }

      ind.traj38 <- c()

      if(check.cdr3.mode==1){
        #Correct TRAJ38 in human (e.g., issue with some 10X data)
        if(species=="HomoSapiens" & ch=="TRA"){

          ind.traj38 <- which((input[ind.species,J]=="TRAJ38" | input[ind.species,J]=="TRAJ38*01") & str_sub(input[ind.species,cdr3],start=-2)=="LI" & nchar(input[ind.species,cdr3]) < Lmax)
          input[ind.species[ind.traj38],cdr3] <- paste(input[ind.species[ind.traj38],cdr3], "W", sep="")

        }

        #Extract the first (start.lg) and last (end.lg) amino acids
        first <- substr(input[ind.species,cdr3], 1, start.lg)
        last <- str_sub(input[ind.species,cdr3], start=-end.lg)

        #print(first)
        #print(last)

        #Find cases incompatible with the reference (allowing matching to any allele)
        nm <- input[ind.species,V]
        rf <- sapply(ref.cdr3.first[[species]][[ch]], function(x){unique(substr(x,1,start.lg))}) #This includes all non-redundant allelic variants

        diff.first <- sapply(1:length(first), function(i){
          diff <- F
          if(!is.na(nm[i]) & !is.na(first[i])){
            diff <- T
            #Check if 'first' matches one of the possible allele of the V gene
            for(st in rf[[nm[i]]]){
              if(nchar(first[i])==nchar(st)){
                if(first[i] == st ){
                  diff <- F
                }
              } else if (nchar(first[i]) > nchar(st)) {
                #This is the special case were 'first' is actually longer than some rf[[nm[i]]]
                #Typically, this is the case when using large values for start.lg
                if(substr(first[i],1,nchar(st)) == st){
                  diff <- F
                }
              }
            }
          }
          return(diff)
        })
        ind.first <- (1:length(first))[diff.first]

        #Find cases incompatible with the reference (allowing matching to any allele)
        nm <- input[ind.species,J]
        rf <- sapply(ref.cdr3.last[[species]][[ch]], function(x){unique(str_sub(x,start=-end.lg))})

        diff.last <- sapply(1:length(last), function(i){
          diff <- F
          if(!is.na(nm[i]) & !is.na(last[i])){
            diff <- T
            #Check if 'last' matches one of the possible allele of the V gene
            for(st in rf[[nm[i]]]){
              if(nchar(last[i])==nchar(st)){
                if(last[i] == st ){
                  diff <- F
                }
              } else if (nchar(last[i]) > nchar(st)) {
                #This is the special case were 'last' is actually longer than some rf[[nm[i]]]
                #Typically, this is the case when using large values for start.lg
                if(str_sub(last[i],start=-nchar(st)) == st){
                  diff <- F
                }
              }
            }
          }
          return(diff)
        })
        ind.last <- (1:length(last))[diff.last]

      } else if (check.cdr3.mode == 2) {


        ## helpers (define once, reused) --------------------------------
        len_col <- function(n) paste0("len_", max(7L, min(22L, as.integer(n))))

        parse_cell <- function(x) {
          if (is.null(x) || is.na(x)) return(character(0))
          x <- as.character(x)
          x <- gsub("\\[|\\]|^c\\(|\\)$|\"", "", x)  # strip brackets / c()
          x <- gsub("^'+|'+$", "", x)                # outer quotes
          parts <- unlist(strsplit(x, "\\s*[,'|\\s]\\s*"))
          parts <- parts[nzchar(parts)]
          unique(parts)
        }

        .first_mismatch_from_start <- function(s, alts) {
          if (length(alts) == 0L || is.na(s)) return(NA_integer_)
          best <- 0L
          for (a in alts) {
            if (!nzchar(a)) next
            m <- min(nchar(s), nchar(a)); k <- 0L
            while (k < m &&
                   substr(s, k + 1L, k + 1L) ==
                   substr(a, k + 1L, k + 1L)) {
              k <- k + 1L
            }
            if (k > best) best <- k
          }
          best + 1L
        }

        .last_mismatch_from_end <- function(s, alts) {
          if (length(alts) == 0L || is.na(s)) return(NA_integer_)
          best <- 0L
          for (a in alts) {
            if (!nzchar(a)) next
            m <- min(nchar(s), nchar(a)); k <- 0L
            while (k < m &&
                   substr(s, nchar(s) - k, nchar(s) - k) ==
                   substr(a, nchar(a) - k, nchar(a) - k)) {
              k <- k + 1L
            }
            if (k > best) best <- k
          }
          -(best + 1L)
        }

        ## lookup tables for this species/chain -------------------------
        TV <- T_V[[species]][[ch]]
        TJ <- T_J[[species]][[ch]]

        # ensure rownames are V/J gene names if needed
        if (is.null(rownames(TV)) && "gene" %in% names(TV)) rownames(TV) <- TV$gene
        if (is.null(rownames(TJ)) && "gene" %in% names(TJ)) rownames(TJ) <- TJ$gene

        ## data for this species subset ---------------------------------
        seqs   <- input[ind.species, cdr3]
        vgenes <- input[ind.species, V]
        jgenes <- input[ind.species, J]
        Ls     <- nchar(seqs)
        cols   <- vapply(Ls, len_col, character(1))

        n <- length(seqs)

        ## precreate output columns on full input -----------------------
        col_consistency <- paste0(ch, "_pred_consistent")
        col_comments    <- paste0("Comments_", ch)

        if (!col_consistency %in% names(input)) {
          input[[col_consistency]] <- NA
        }
        if (!col_comments %in% names(input)) {
          input[[col_comments]] <- vector("list", nrow(input))
        }

        # indices where checking V side makes sense
        useV <- !is.na(seqs) & !is.na(vgenes) &
          vgenes %in% rownames(TV) & cols %in% colnames(TV)

        keysV <- paste(vgenes[useV], cols[useV], sep = "|")
        ukeysV <- unique(keysV)
        cacheV <- setNames(vector("list", length(ukeysV)), ukeysV)
        for (k in ukeysV) {
          parts <- strsplit(k, "\\|", fixed = FALSE)[[1]]
          g  <- parts[1]
          lc <- parts[2]
          cacheV[[k]] <- parse_cell(TV[g, lc])
        }
        allowedV_all <- vector("list", n)
        allowedV_all[useV] <- cacheV[keysV]

        # indices where checking J side makes sense
        useJ <- !is.na(seqs) & !is.na(jgenes) &
          jgenes %in% rownames(TJ) & cols %in% colnames(TJ)

        keysJ <- paste(jgenes[useJ], cols[useJ], sep = "|")
        ukeysJ <- unique(keysJ)
        cacheJ <- setNames(vector("list", length(ukeysJ)), ukeysJ)
        for (k in ukeysJ) {
          parts <- strsplit(k, "\\|", fixed = FALSE)[[1]]
          g  <- parts[1]
          lc <- parts[2]
          cacheJ[[k]] <- parse_cell(TJ[g, lc])
        }
        allowedJ_all <- vector("list", n)
        allowedJ_all[useJ] <- cacheJ[keysJ]


        diff.first <- rep(NA, n)
        diff.last  <- rep(NA, n)

        comments_local <- vector("list", n)  # comments per row in species subset

        for (i in seq_len(n)) {

          s  <- seqs[i]
          allowedV <- allowedV_all[[i]]
          allowedJ <- allowedJ_all[[i]]

          v_comment <- ""
          j_comment <- ""

          ## V-side decision
          if (length(allowedV) == 0L || is.na(s)) {
            # no check -> keep NA in diff.first
            diff.first[i] <- NA
          } else {
            ok <- any(vapply(
              allowedV,
              function(a) {
                nzchar(a) &&
                  nchar(s) >= nchar(a) &&
                  substr(s, 1L, nchar(a)) == a
              },
              logical(1)
            ))
            if (ok) {
              diff.first[i] <- FALSE
            } else {
              pos <- .first_mismatch_from_start(s, allowedV)
              v_comment <- paste("V: Inconsistency at position", pos)
              diff.first[i] <- TRUE
            }
          }

          ## J-side decision
          if (length(allowedJ) == 0L || is.na(s)) {
            diff.last[i] <- NA
          } else {
            ok <- any(vapply(
              allowedJ,
              function(a) {
                nzchar(a) &&
                  nchar(s) >= nchar(a) &&
                  substr(s, nchar(s) - nchar(a) + 1L, nchar(s)) == a
              },
              logical(1)
            ))
            if (ok) {
              diff.last[i] <- FALSE
            } else {
              pos <- .last_mismatch_from_end(s, allowedJ)
              j_comment <- paste("J: Inconsistency at position", pos)
              diff.last[i] <- TRUE
            }
          }

          # store comments for this row (if any)
          if (nzchar(v_comment) || nzchar(j_comment)) {
            comm <- c(v_comment, j_comment)
            comm <- comm[nzchar(comm)]                      # drop empty ones
            comments_local[[i]] <- paste(comm, collapse = "; ")
          }
        }

        ## ------------------------------------------------------------
        ## 3) Build consistent / inconsistent with NA rules
        ## ------------------------------------------------------------

        both_na  <- is.na(diff.first) & is.na(diff.last)
        any_true <- (diff.first %in% TRUE) | (diff.last %in% TRUE)

        inconsistent <- ifelse(both_na, NA, any_true)
        consistent   <- ifelse(is.na(inconsistent), NA, !inconsistent)

        input[ind.species, col_consistency] <- consistent

        ## indices used later in verbose output (ignore NAs)
        ind.first <- which(!is.na(diff.first) & diff.first)
        ind.last  <- which(!is.na(diff.last)  & diff.last)

        ## write comments back in one go for this species subset
        has_comm <- !vapply(comments_local, is.null, logical(1))
        if (any(has_comm)) {
          input[[col_comments]][ind.species[has_comm]] <- comments_local[has_comm]
        }
      }
    }
  }

  return(input)
}


correct.VJnames <- function(input, segment.list=c("TRAV","TRAJ","TRBV","TRBJ"), species.default="HomoSapiens",
                            use.allele=F, seq.protocol="Default", verbose=1){

  if("species" %in% colnames(input)){
    species.list <- unique(input[,"species"])
  } else {
    species.list <- c(species.default)
  }


  for(species in species.list){
    for(s in segment.list){
      if(use.allele){
        name.list <- gene.allele.list[[seq.protocol]][[species]][substr(gene.allele.list[[seq.protocol]][[species]],1,4)==s]
      } else {
        name.list <- gene.list[[seq.protocol]][[species]][substr(gene.list[[seq.protocol]][[species]],1,4)==s]
      }
      if("species" %in% colnames(input)){
        ind <- which(!(input[,s] %in% name.list) & input[,"species"]==species & !is.na(input[,s]))
      } else {
        ind <- which(!(input[,s] %in% name.list) & !is.na(input[,s]))
      }

      if(length(ind)>0){

        nm <- strsplit(input[ind,s], split="*", fixed=T)
        gene <- unlist(lapply(nm, function(x){x[1]}))
        allele <- unlist(lapply(nm, function(x){x[2]}))

        ga <- unlist(lapply(1:length(gene), function(x){ clean.name.allele(gene=gene[x], allele=allele[x], species=species,
                                                                           use.allele=use.allele, seq.protocol=seq.protocol)}))

        #Set to NA cases where the gene names comes from another segment
        #This is because the clean.name.allele does not check if the gene is in the right column (e.g., TRAV12-2 in TRAJ column)
        ga[substr(ga,1,4) != s] <- NA

        if(verbose>0){
          i <- which(input[ind,s] != ga & is.na(ga)==F)
          if(length(i)>0){

            m.cor <- data.frame(original.name = input[ind[i],s], corrected.name = ga[i],row.names = NULL)
            m.cor <- m.cor[!duplicated(m.cor),]

            entry <- ifelse(length(i)==1,"entry","entries")
            nm <- ifelse(dim(m.cor)[1]==1,"name","names")
            verb <- ifelse(dim(m.cor)[1]==1,"was","were")
            print(paste("*** ",dim(m.cor)[1]," ",s," ",nm," in ",length(i)," ",entry," ",verb," corrected in ",species, "***",sep=""))
            if(verbose==1 | verbose==2){
              print("Use verbose=3 to see them")
            }
            if(verbose==3){
              print(m.cor)
            }

            cat("\n")
          }

          #Check the cases where the segment was not NA, but was put to NA (i.e., mapping of gene name failed)
          i <- which(!is.na(input[ind,s]) & input[ind,s] != "" & is.na(ga)==T)

          if(length(i)>0){
            v <- unique(input[ind[i],s])
            v <- v[!is.na(v)]
            print(paste("*** ",length(v), " ", s, " names in ",length(i)," entries could not be corrected in ",species," ***", sep=""))
            if(verbose==1){
              n <- min(10,length(v))
              if(length(v)>n){
                print("Examples  (use verbose >= 2 to see them all):")
              }
              print(v[1:n])
            }
            if(verbose>1){
              print(v)
            }
            cat("\n")
          }
        }

        input[ind,s] <- ga
      }
    }
  }
  return(input)
}

verify.chain <- function(input, chain){

  # Check cases where people leave the chain="AB",
  # but actually provide single chain data, including with ambiguous colnames like V,J,CDR3_seq

  format <- determine.format(input, verbose=0)
  # Only check for the format not based on clone.id (cases with formats based on clone.id will always be treated as alpha+beta, even if one chain is empty)
  if(!format %in% names(clone.format.col)){
    if(chain=="AB"){

      map.A <- names(mapping.colnames[["AB"]][which(mapping.colnames[["AB"]] %in% c("TRAV", "TRAJ", "cdr3_TRA"))])
      map.B <- names(mapping.colnames[["AB"]][which(mapping.colnames[["AB"]] %in% c("TRBV", "TRBJ", "cdr3_TRB"))])
      inter.A <- intersect(colnames(input), c("TRAV", "TRAJ", "cdr3_TRA", map.A))
      inter.B <- intersect(colnames(input), c("TRBV", "TRBJ", "cdr3_TRB", map.B))

      if(length(inter.A)==0 | length(inter.B)==0){
        if(length(inter.A)==0 & length(inter.B)>0){
          print("Missing columns for alpha chain, only beta chain will be considered")
          chain <- "B"
        } else if(length(inter.A)>0 & length(inter.B)==0){
          print("Missing columns for beta chain, only alpha chain will be considered")
          chain <- "A"
        } else if(length(inter.A)==0 & length(inter.B)==0){
          # Check if the a single chain can be inferred.
          # This is the case for instance when providing data in AIRR format without specifying the chain
          map.unk <- names(mapping.colnames[["A"]][which(mapping.colnames[["A"]] == "TRAV")])
          if(length(intersect(map.unk, colnames(input)))==1){
            p <- intersect(map.unk, colnames(input))[1]
            if("TRA" %in% unique(substr(input[,p],1,3)) & ! "TRB" %in% unique(substr(input[,p],1,3))){
              print("Missing data for beta chain, only alpha chain will be considered")
              chain <- "A"
            } else if("TRB" %in% unique(substr(input[,p],1,3)) & ! "TRA" %in% unique(substr(input[,p],1,3))){
              print("Missing data for alpha chain, only beta chain will be considered")
              chain <- "B"
            } else {
              chain <- ""
            }
          }
        }
      }
    } else if(chain=="A"){

      map.A <- names(mapping.colnames[["A"]][which(mapping.colnames[["A"]] %in% c("TRAV", "TRAJ", "cdr3_TRA"))])
      inter.A <- intersect(colnames(input), c("TRAV", "TRAJ", "cdr3_TRA", map.A))
      if(length(inter.A)==0){
        chain <- ""
      }

    } else if(chain=="B"){

      map.B <- names(mapping.colnames[["B"]][which(mapping.colnames[["B"]] %in% c("TRBV", "TRBJ", "cdr3_TRB"))])
      inter.B <- intersect(colnames(input), c("TRBV", "TRBJ", "cdr3_TRB", map.B))
      if(length(inter.B)==0){
        chain <- ""
      }

    }
  }
  return(chain)
}
#Determine the format
determine.format <- function(input, verbose=1){

  col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")
  col.A <- col[1:3]
  col.B <- col[4:6]

  format.list <- names(clone.format.col)
  format <- "custom"
  if(length(intersect(colnames(input), col.A))==3 | length(intersect(colnames(input), col.B))==3){
    format <- "MixTCRclean"
  } else {
    for(f in format.list){
      #Check that the three field for specific formats are present, and none of the standard colnames in MixTCRclean
      if(length(intersect(colnames(input), clone.format.col[[f]]))==3 & length(intersect(colnames(input),col))==0){
        format <- f
        if(verbose==1){
          print(paste("Inferred format:",f))
        }
        break
      }
    }
  }

  return(format)
}

stack_clones <- function(input, format){

  if(format %in% names(clone.format.col)){

    col <- clone.format.col[[format]]
    other.col <- setdiff(colnames(input), col)

    input.f <- apply(input,1,function(x){
      if(substr(x[col[1]],1,3)=="TRA"||substr(x[col[1]],1,4)=="TCRA"){
        v <- c(x[col],NA,NA,NA,x[other.col])
      } else if(substr(x[col[1]],1,3)=="TRB"||substr(x[col[1]],1,4)=="TCRB"){
        v <- c(NA,NA,NA,x[col],x[other.col])
      } else if(format=="VDJdb" & "Gene" %in% names(x)){
        if(x["Gene"]=="TRA"){
          v <- c(x[col],NA,NA,NA,x[other.col])
        } else if(x["Gene"]=="TRB"){
          v <- c(NA,NA,NA,x[col],x[other.col])
        }
      } else if("chain" %in% names(x)){
        if(x["chain"]=="TRA"){
          v <- c(x[col],NA,NA,NA,x[other.col])
        } else if(x["chain"]=="TRB"){
          v <- c(NA,NA,NA,x[col],x[other.col])
        }
      } else {
        #In this case, inferring the chain failed
        v <- c(NA,NA,NA,NA,NA,NA,x[other.col])
      }
      return(v)
    })
    input.f <- t(input.f)
    colnames(input.f) <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB",other.col)
    input.f <- as.data.frame(input.f)
    if(format=="VDJdb" & "complex.id" %in% colnames(input.f)){
      input.f[,"complex.id"] <- as.numeric(input.f[,"complex.id"])
    }
  } else {
    input.f <- input
  }
  return(input.f)
}

clean.name.allele <- function(gene=gene, allele=allele, species="HomoSapiens", seq.protocol="Default", use.allele=F){

  if(species != "HomoSapiens" & species != "MusMusculus"){
    print("Undefined species: ",species)
  }

  if(!use.allele){
    allele <- ""
  }

  if(is.na(gene) | gene==""){  #This is not needed in MixTCRclean, but can be useful in other cases
    ga <- NA
  } else {

    #Check if the gene needs to be corrected
    if(gene %in% gene.list[[seq.protocol]][[species]] == F){

      #Do a few automatic corrections
      gene <- gsub("TCR","TR",gene)
      gene <- gsub("–","-", gene)
      gene <- gsub("-0","-", gene)
      gene <- sub("-$", "", gene)
      if(species=="HomoSapiens"){ gene <- gsub("hTR","TR",gene) }
      if(species=="MusMusculus"){ gene <- gsub("mTR","TR",gene) }
      st0 <- substr(gene,1,2)
      st <- substr(gene,3,4)
      if(st0=="TR" & st %in% c("AV","AJ","BV","BJ")){
        gene <- gsub(paste("TR",st,"0",sep=""),paste("TR",st,sep=""),gene)
      }
      if(substr(gene,1,4)=="TRAJ" & grepl("-",gene)){
        gene <- unlist(strsplit(gene,"-"))[1]
      }

      if(gene %in% gene.list[[seq.protocol]][[species]] == F ){
        #The gene is still not ok, try using our manual dictionary
        if(!is.na(map[[species]][gene])){
          gene <- map[[species]][gene] #The gene was wrong, but can be corrected

        } else {
          gene <- NA #The gene was wrong and could not be corrected
        }
      }
    }

    if(use.allele){
      if(!is.na(gene)){
        #The gene was ok or could be corrected
        v <- paste(gene,allele,sep="*")
        if(v %in% gene.allele.list[[seq.protocol]][[species]] == T){
          ga <- v #the gene*allele is ok
        } else {
          ga <- paste(gene,allele.default[[seq.protocol]][[species]][gene],sep="*") #If not, Use the gene*default.allele
        }
      } else {
        #The gene could not be corrected
        ga <- NA
      }
    } else {
      ga <- gene
    }
  }
  return(ga)
}

#' @export
merge_mouse_TRAV <- function(input){

  # This has to be run after alleles have been removed and genes have been corrected
  # If the "species" field is present, it takes only "MusMusculus" entries
  # If not, it assumes all entries are "MusMusculus"
  # It also assumes that alleles have been removed

  if("TRAV" %in% colnames(input)){

    if("species" %in% colnames(input)){
      ind <- which(input[,"species"]=="MusMusculus")
    } else {
      ind <- 1:dim(input)[1]
    }

    v.cor <- as.character(unlist(lapply(input[["TRAV"]][ind], function(y){  # WARNING: I don't understand why not using es[ind,"TRAV"]
      if (y %in% names(merge.mouse.TRAV)){
        y <- merge.mouse.TRAV[y]
      }
      return(y)
    })))

    input[ind,"TRAV"] <- v.cor
  }
  return(input)

}

#Build actual clones if the data are provided in a format based on clone IDs.
merge_clones <- function(input, format){

  col <- c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB")

  if(length(intersect(colnames(input), clone.id))==1){

    TCR.col <- clone.format.col[[format]]
    Vn <- TCR.col[1]
    barcode.label <- intersect(colnames(input), clone.id)

    #other.col <- setdiff(colnames(input), TCR.col) #keep all columns (can be complex if different values are used with the same barcode)
    other.col <- intersect(colnames(input),c(barcode.label,"model","species")) #keep only the barcode, the model and the species (if present)
    input.f <- as.data.frame(matrix(nrow=0, ncol=6+length(other.col)))

    if("model" %in% colnames(input)){
      barcode.list <- paste(input[,barcode.label], input[,"model"]) #Include the possibility of having redundant barcodes for different models.
    } else {
      barcode.list <- input[,barcode.label]
    }

    input.t <- input[order(barcode.list),]
    barcode.list <- sort(barcode.list)

    #Go through all clone IDs
    i <- 1
    while(i<=dim(input.t)[1]){
      seq.A <- list()
      seq.B <- list()
      ct.A <- 0
      ct.B <- 0
      #This is the special case of the '0' complex.id for VDJdb which contains all single chain entries
      if(input.t[i,barcode.label]==0 & format=="VDJdb"){
        ind.A <- which(input.t[,barcode.label]==0 & substr(input.t$V,1,3)=="TRA")
        ind.B <- which(input.t[,barcode.label]==0 & substr(input.t$V,1,3)=="TRB")
        if(length(ind.A)>0){
          input.f <- rbind(input.f,cbind(input.t[ind.A,TCR.col],NA,NA,NA,input.t[ind.A,other.col]))
        }
        if(length(ind.B)>0){
          input.f <- rbind(input.f,cbind(NA,NA,NA,input.t[ind.B,TCR.col],input.t[ind.B,other.col]))
        }
        i <- i+length(ind.A)+length(ind.B)

      } else {
        for(j in 0:(dim(input.t)[1]-i)){
          if(barcode.list[i+j]==barcode.list[i]){
            if(substr(input.t[i+j,Vn],1,3)=="TRA"){
              ct.A <- ct.A+1
              seq.A[[ct.A]] <- as.character(input.t[i+j,TCR.col])
            } else if(substr(input.t[i+j,Vn],1,3)=="TRB"){
              ct.B <- ct.B+1
              seq.B[[ct.B]] <- as.character(input.t[i+j,TCR.col])
            }
          } else {
            break
          }
        }


        if(length(seq.A)>0 & length(seq.B)>0){
          for(sA in seq.A){
            for(sB in seq.B){
              seq.all <- c(sA,sB,unlist(input.t[i,other.col]))
              input.f <- rbind(input.f, seq.all)
            }
          }
        } else if(length(seq.A)>0 & length(seq.B)==0){
          for(sA in seq.A){
            seq.all <- c(sA,c(NA,NA,NA),unlist(input.t[i,other.col]))
            input.f <- rbind(input.f, seq.all)
          }
        } else if(length(seq.A)==0 & length(seq.B)>0){
          for(sB in seq.B){
            seq.all <- c(c(NA,NA,NA),sB,unlist(input.t[i,other.col]))
            input.f <- rbind(input.f, seq.all)
          }
        }
        i <- i+ct.A+ct.B
      }
    }
    colnames(input.f) <- c(col, other.col)
    if(format=="VDJdb" & "complex.id" %in% colnames(input.f)){
      input.f[,"complex.id"] <- as.numeric(input.f[,"complex.id"])
      print(input.f[,"complex.id"])
    }
    return(input.f)

  } else {
    print(paste("WARNING: Unable to reconstruct clones. Clone_id should be exactly one element in",paste(clone.id, collapse=" / ")))
    print("  Try using build.clones = F")
    return(input)
  }


}



