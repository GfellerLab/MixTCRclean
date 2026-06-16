
#####
# Define some function
#####


#' @export
clean_input.MixTCRclean <- function(input, use.allele=F, correct.gene.names=T, use.mouse.strain=F,
                                    chain="AB", species.default="HomoSapiens", check.cdr3.mode=2,
                                    keep.incomplete.chain=T, start.lg=1, end.lg=2, seq.protocol="IMGT",
                                    merge.ambiguous=T, verbose=1){

  ####
  # Clean the input by removing CDR3 with weird characters, longer than Lmax or shorter than Lmin
  # Correct VJ genes based on our dictionary
  # species.default is only used if input does not contain the "species" column
  # merge.ambiguous should be set to FALSE ONLY if clean_input is used outside of MixTCRclean
  # this will prevent mapping TRBV6-2 to TRBV6-2/6-3
  ####

  if(seq.protocol=="SEQTR"){
    merge.TRBV6_2_3 <- T
    merge.TRBV12_3_4 <- T
  }
  if(seq.protocol=="Default"){
    merge.TRBV6_2_3 <- T
    merge.TRBV12_3_4 <- F
  }
  if(seq.protocol=="IMGT"){
    merge.TRBV6_2_3 <- F
    merge.TRBV12_3_4 <- F
  }

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
    ind <- which( nc < MixTCRviz::Lmin | nc >  MixTCRviz::Lmax | grepl('[^ACDEFGHIKLMNPQRSTVWY]', input[,cdr3]) == T)
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

  if(merge.TRBV6_2_3){
    #Do a manual correction for TRBV6-2 and TRBV6-3 -> TRBV6-2/6-3
    if(chain=="B" | chain=="AB"){

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
  if(merge.TRBV12_3_4){
    #Do a manual correction for TRBV12-3 and TRBV12-4 -> TRBV12-3/12-4
    #This is always done if seq.protocol is SEQTR
    if(chain=="B" | chain=="AB"){

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
    }
  }

  # Correct gene names
  # If alleles, it will correct the gene name, and keep the allele. If the allele cannot be found, it will remove it
  # If genes, it will correct the gene name
  # If gene name cannot be corrected, it gives NA

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
    input <- check_cdr3.MixTCRclean(input=input, chain=chain, species.default=species.default, check.cdr3.mode=check.cdr3.mode,
                        start.lg=start.lg, end.lg=end.lg, verbose=verbose)
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
check_cdr3.MixTCRclean <- function(input, chain="AB", species.default="HomoSapiens", check.cdr3.mode=2, start.lg=1, end.lg=2, verbose=1){

  # Clean the CDR3 based on the V and J usage.
  # This should be applied after correcting the gene names, and adding the species if needed
  # species.default is only used if es.all does not contain the "species" column
  # If the allele is given in the gene name, the allele will be used.
  chain.list <- paste0("TR",strsplit(chain,split="")[[1]])
  chain.small <- tolower(strsplit(chain,split="")[[1]])
  names(chain.small) <- chain.list

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
          ind.traj38 <- which((input[ind.species,J]=="TRAJ38" | input[ind.species,J]=="TRAJ38*01") & str_sub(input[ind.species,cdr3],start=-2)=="LI" & nchar(input[ind.species,cdr3]) <  MixTCRviz::Lmax)
          input[ind.species[ind.traj38],cdr3] <- paste(input[ind.species[ind.traj38],cdr3], "W", sep="")
        }

        ## --- Setup output columns (same as mode 2) ----------------------
        col_consistency <- paste0(ch, "_pred_consistent")
        col_comments    <- paste0("Comments_", ch)

        if (!col_consistency %in% names(input)) {
          input[[col_consistency]] <- NA
        }
        if (!col_comments %in% names(input)) {
          input[[col_comments]] <- vector("list", nrow(input))
        }

        ## --- Helpers to locate the first mismatch -----------------------
        .first_mismatch_from_start <- function(s, alts) {
          best <- 0L
          for (a in alts) {
            if (!nzchar(a)) next
            m <- min(nchar(s), nchar(a)); k <- 0L
            while (k < m && substr(s, k+1L, k+1L) == substr(a, k+1L, k+1L)) k <- k + 1L
            if (k > best) best <- k
          }
          best + 1L
        }
        .last_mismatch_from_end <- function(s, alts) {
          best <- 0L
          for (a in alts) {
            if (!nzchar(a)) next
            m <- min(nchar(s), nchar(a)); k <- 0L
            while (k < m &&
                   substr(s, nchar(s)-k, nchar(s)-k) ==
                   substr(a, nchar(a)-k, nchar(a)-k)) k <- k + 1L
            if (k > best) best <- k
          }
          -(best + 1L)
        }

        ## --- Extract the first (start.lg) and last (end.lg) amino acids ---
        first <- substr(input[ind.species,cdr3], 1, start.lg)
        last  <- str_sub(input[ind.species,cdr3], start=-end.lg)
        n     <- length(first)

        ## --- V-side check -----------------------------------------------
        nm.v <- input[ind.species,V]
        rf.v <- sapply(MixTCRviz::ref.cdr3.first[[species]][[ch]],
                       function(x){ unique(substr(x,1,start.lg)) })

        diff.first <- sapply(1:n, function(i){
          diff <- F
          if(!is.na(nm.v[i]) & !is.na(first[i])){
            diff <- T
            for(st in rf.v[[nm.v[i]]]){
              if(nchar(first[i])==nchar(st)){
                if(first[i] == st){ diff <- F }
              } else if (nchar(first[i]) > nchar(st)) {
                if(substr(first[i],1,nchar(st)) == st){ diff <- F }
              }
            }
          }
          return(diff)
        })
        ind.first <- (1:n)[diff.first]

        ## --- J-side check -----------------------------------------------
        nm.j <- input[ind.species,J]
        rf.j <- sapply(MixTCRviz::ref.cdr3.last[[species]][[ch]],
                       function(x){ unique(str_sub(x,start=-end.lg)) })

        diff.last <- sapply(1:n, function(i){
          diff <- F
          if(!is.na(nm.j[i]) & !is.na(last[i])){
            diff <- T
            for(st in rf.j[[nm.j[i]]]){
              if(nchar(last[i])==nchar(st)){
                if(last[i] == st){ diff <- F }
              } else if (nchar(last[i]) > nchar(st)) {
                if(str_sub(last[i],start=-nchar(st)) == st){ diff <- F }
              }
            }
          }
          return(diff)
        })
        ind.last <- (1:n)[diff.last]

        ## --- Build comments and consistency column ----------------------
        seqs <- input[ind.species, cdr3]
        comments_local <- vector("list", n)

        for (i in seq_len(n)) {
          v_comment <- ""
          j_comment <- ""

          if (diff.first[i]) {
            pos <- .first_mismatch_from_start(seqs[i], rf.v[[nm.v[i]]])
            v_comment <- paste("V: Inconsistency at position", pos)
          }
          if (diff.last[i]) {
            pos <- .last_mismatch_from_end(seqs[i], rf.j[[nm.j[i]]])
            j_comment <- paste("J: Inconsistency at position", pos)
          }
          if (nzchar(v_comment) || nzchar(j_comment)) {
            comm <- c(v_comment, j_comment)
            comm <- comm[nzchar(comm)]
            comments_local[[i]] <- paste(comm, collapse = "; ")
          }
        }

        ## consistency column: TRUE if neither side flagged a mismatch
        input[ind.species, col_consistency] <- !(diff.first | diff.last)

        ## write comments back in one go
        has_comm <- !vapply(comments_local, is.null, logical(1))
        if (any(has_comm)) {
          input[[col_comments]][ind.species[has_comm]] <- comments_local[has_comm]
        }
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

        # both_na  <- is.na(diff.first) & is.na(diff.last)
        # any_true <- (diff.first %in% TRUE) | (diff.last %in% TRUE)
        #
        # inconsistent <- ifelse(both_na, NA, any_true)
        # consistent   <- ifelse(is.na(inconsistent), NA, !inconsistent)

        any_true   <- (diff.first %in% TRUE) | (diff.last %in% TRUE)
        consistent <- ifelse(is.na(seqs), NA, !any_true)

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
  for (ch in chain.list) {
    col_comments <- paste0("Comments_", ch)
    if (col_comments %in% names(input)) {
      input[[col_comments]] <- vapply(input[[col_comments]], function(x) {
        if (is.null(x) || (length(x) == 1L && is.na(x))) "" else as.character(x)
      }, character(1))
    }
  }
  return(input)
}


