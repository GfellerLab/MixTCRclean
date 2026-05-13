# Description ---------------------------------------------------------------
# In this file we create and export some variables used by the package
# The script should be run from MixTCRviz folder

library(ggplot2)
library(stringr)

# Mapping of the gene names -----------------------------------------------
# These are only access via MixTCRviz:::
# segment.list <- c("TRAV", "TRBV", "TRAJ", "TRBJ")
# species.list <- c("HomoSapiens", "MusMusculus")
#
# # Take the list of alleles
#
#
# # gene.allele.list <- list(); gene.allele.list[["IMGT"]] <- list()
# # gene.list <- list(); gene.list[["IMGT"]] <- list()
# # gene.type <- list(); gene.type[["IMGT"]] <- list()
# # allele.default <- list(); allele.default[["IMGT"]] <- list()
# #
# # for(tsp in species.list){
# #   gene.allele.list$IMGT[[tsp]] <- c()
# #   gene.list$IMGT[[tsp]] <- c()
# #   gene.type$IMGT[[tsp]] <- c()
# #   for(s in segment.list){
# #     df <- read.csv(file=paste("../MixTCRviz/data_raw/CDR123/",tsp,"/",s,"_allele.csv", sep=""), row.names = 1)
# #     gene.allele.list$IMGT[[tsp]] <- c(gene.allele.list$IMGT[[tsp]], rownames(df))
# #     gt <- df[,"gene_type"]
# #     names(gt) <- rownames(df)
# #     gene.type$IMGT[[tsp]] <- c(gene.type$IMGT[[tsp]], gt)
# #
# #     df <- read.csv(file=paste("../MixTCRviz/data_raw/CDR123/",tsp,"/",s,".csv", sep=""), row.names = 1)
# #     gene.list$IMGT[[tsp]] <- c(gene.list$IMGT[[tsp]], rownames(df))
# #
# #     gt <- df[,"gene_type"]
# #     names(gt) <- rownames(df)
# #     gene.type$IMGT[[tsp]] <- c(gene.type$IMGT[[tsp]], gt)
# #
# #     default <- df[,"default_allele"]
# #     default <- gsub('\\*','',default)
# #     names(default) <- rownames(df)
# #     allele.default$IMGT[[tsp]] <- c(allele.default$IMGT[[tsp]], default) #This is no longer used.
# #   }
# # }
# #
#
#
# #For mouse TRAV, take the merging between different strains
# m <- read.csv(paste("data_raw/CDR123/MusMusculus/TRAV_merged.csv", sep=""))
# merge.mouse.TRAV <- m[,2]
# names(merge.mouse.TRAV) <- m[,1]
#
# mp.TRAV <- read.csv("data_raw/TidyVJ/mapping_TRAV.csv", header=F, skip=1)
# mp.TRBV <- read.csv("data_raw/TidyVJ/mapping_TRBV.csv", header=F, skip=1)
# mp.TRAJ <- read.csv("data_raw/TidyVJ/mapping_TRAJ.csv", header=F, skip=1)
# mp.TRBJ <- read.csv("data_raw/TidyVJ/mapping_TRBJ.csv", header=F, skip=1)
# mp <- rbind(mp.TRAV,mp.TRBV,mp.TRAJ,mp.TRBJ)
# # Keep only the entries that could be mapped
# mp <- mp[mp[,5]=="yes",]
#
# #Build the mapping separately for human and mouse
# map <- list()
# pos.h <- which(mp[,4]=="HomoSapiens")
# map[["HomoSapiens"]] <- mp[pos.h,2]
# names(map[["HomoSapiens"]]) <- mp[pos.h,1]
# pos.m <- which(mp[,4]=="MusMusculus")
# map[["MusMusculus"]] <- mp[pos.m,2]
# names(map[["MusMusculus"]]) <- mp[pos.m,1]
#
# #For SEQTR, we merge TRBV12-3 and TRBV12-4 into TRBV12-3/12-4
#
# gene.allele.list[["Default"]] <- gene.allele.list$IMGT
# gene.list[["Default"]] <-  gene.list$IMGT
# gene.type[["Default"]] <-  gene.type$IMGT
# allele.default[["Default"]] <-  allele.default$IMGT
#
# gene.Default.exclude <- c("TRBV6-2*01", "TRBV6-3*01", "TRBV6-2", "TRBV6-3")
#
# ind <- which(! gene.allele.list$Default$HomoSapiens %in% gene.Default.exclude)
# gene.allele.list$Default$HomoSapiens <- sort(c(gene.allele.list$Default$HomoSapiens[ind], "TRBV6-2/6-3*01"))
# ind <- which(! gene.list$Default$HomoSapiens %in% gene.Default.exclude)
# gene.list$Default$HomoSapiens <- sort(c(gene.list$Default$HomoSapiens[ind], "TRBV6-2/6-3"))
# ind <- which(! names(gene.type$Default$HomoSapiens) %in% gene.Default.exclude)
# gene.type$Default$HomoSapiens <- setNames(c(gene.type$Default$HomoSapiens[ind], "F"),c(names(gene.type$Default$HomoSapiens[ind]), "TRBV6-2/6-3"))
# ind <- which(! names(allele.default$Default$HomoSapiens) %in% gene.Default.exclude)
# allele.default$Default$HomoSapiens <- setNames(c(allele.default$Default$HomoSapiens[ind], "01"),c(names(allele.default$Default$HomoSapiens[ind]), "TRBV6-2/6-3"))
#
#
# gene.allele.list[["SEQTR"]] <- gene.allele.list$Default
# gene.list[["SEQTR"]] <-  gene.list$Default
# gene.type[["SEQTR"]] <-  gene.type$Default
# allele.default[["SEQTR"]] <-  allele.default$Default
#
# gene.SEQTR.exclude <- c("TRBV12-3*01", "TRBV12-4*01", "TRBV12-4*02", "TRBV12-3", "TRBV12-4")
#
# ind <- which(! gene.allele.list$SEQTR$HomoSapiens %in% gene.SEQTR.exclude)
# gene.allele.list$SEQTR$HomoSapiens <- sort(c(gene.allele.list$SEQTR$HomoSapiens[ind], "TRBV12-3/12-4*01"))
# ind <- which(! gene.list$SEQTR$HomoSapiens %in% gene.SEQTR.exclude)
# gene.list$SEQTR$HomoSapiens <- sort(c(gene.list$SEQTR$HomoSapiens[ind], "TRBV12-3/12-4"))
# ind <- which(! names(gene.type$SEQTR$HomoSapiens) %in% gene.SEQTR.exclude)
# gene.type$SEQTR$HomoSapiens <- setNames(c(gene.type$SEQTR$HomoSapiens[ind], "F"),c(names(gene.type$SEQTR$HomoSapiens[ind]), "TRBV12-3/12-4"))
# ind <- which(! names(allele.default$SEQTR$HomoSapiens) %in% gene.SEQTR.exclude)
# allele.default$SEQTR$HomoSapiens <- setNames(c(allele.default$SEQTR$HomoSapiens[ind], "01"),c(names(allele.default$SEQTR$HomoSapiens[ind]), "TRBV12-3/12-4"))
#
#
#
#
# # Defining color/shape maps for TCR genes ----------------------------------
# TCRgene2aes <- list()
# segment.list <- c("TRAV", "TRBV", "TRAJ", "TRBJ")
# palette <- "Set 3"
# # Different color for each gene per segment
# for (tsp in species.list){
#   for (s in segment.list){
#     cGenes <- grep(s, gene.list$IMGT[[tsp]], value=TRUE)
#     cols <- hcl.colors(n=length(cGenes), palette=palette)
#     cGenes <- c(cGenes, "Other")
#     cols <- c(cols, "gray90")
#     # Add a light gray for the "Other" corresponding to the sum of the genes
#     # not showed in the bar plot.
#     TCRgene2aes[[tsp]][[s]]$color2 <- setNames(cols, cGenes)
#   }
# }
#
# # Combining color with shape and outer shape color
# for (tsp in species.list){
#   for (s in segment.list){
#     cGenes <- grep(s, gene.list$IMGT[[tsp]], value=TRUE)
#     genesCode <- sort(unique(gsub("TR(A|B)(V|J)", "", cGenes)))
#
#     if(s == "TRBJ" ){
#       #The number of TRBJ genes is much lower, so we use only two shapes
#       shapes <- c(21,25) # Shapes with outer colors can only take values 21:25
#       n_diffShapes <- length(shapes)
#       g1 <- genesCode
#       g2 <- rep("none", length(g1))
#     } else if ((s=="TRBV" & tsp=="MusMusculus")) {
#       shapes <- c(21,24,25) # Shapes with outer colors can only take values 21:25
#       n_diffShapes <- length(shapes)
#       g1 <- genesCode
#       g2 <- rep("none", length(g1))
#     } else {
#       shapes <- 21:25 # Shapes with outer colors can only take values 21:25
#       n_diffShapes <- length(shapes)
#       g1 <- gsub("-(.*)$", "", genesCode)
#       g2 <- stringr::str_replace(string=genesCode, pattern=g1, replacement="")
#       # g1 is the first part of gene names before the "-" and g2 is the 2nd
#       # part of this name (either "" when there isn't such 2nd part, "-1", "-2", ...).
#       g2 <- gsub("^$", "none", g2)
#     }
#     # Replace empty strings by none to be able to call them by name below.
#     g1u <- sort(unique(g1))
#     g2u <- setdiff(sort(unique(g2)), "none")
#     # The cases with "none" will be considered with same shape/color as first case.
#
#     if (length(unique(g2)) == 1 ){
#       # When there weren't any gene names with "-" in them, we'll instead
#       # combine the color, shape and outer color for all the genes (I'll
#       # consider here 2 different outer colors, and the number of inner color
#       # will be determined based on number needed to distinguish all cases).
#       if(s == "TRBJ" | (s=="TRBV" & tsp=="MusMusculus")){
#         n_grays <- 1
#       } else {
#         #The number of TRBJ genes is lower, so we use only a single outer color
#         n_grays <- 2
#       }
#       g1 <- paste0("l", (rep(1:ceiling(length(g1) / (n_grays*n_diffShapes)),
#         each=n_grays*n_diffShapes)[1:length(g1)]))
#       g2 <- paste0("l", (rep(1 : (n_grays*n_diffShapes), length.out=length(g2))))
#       # Use paste to make sure these are characters to avoid issue as we'll use
#       # these as names below.
#       g1u <- unique(g1)
#       g2u <- unique(g2)
#     }
#
#     cols <- hcl.colors(n=length(g1u), palette=palette)
#     names(cols) <- g1u
#     shapes <- rep(shapes, length.out=length(g2u))
#     cols_out <- gray(seq(0, 0.5, length.out=ceiling(length(g2u) / n_diffShapes)))
#     cols_out <- rep(cols_out, each=n_diffShapes)[1:length(g2u)]
#     shapes <- c(shapes[1], shapes)
#     cols_out <- c(cols_out[1], cols_out)
#     names(shapes) <- names(cols_out) <- c("none", g2u)
#
#     TCRgene2aes[[tsp]][[s]]$color1 <- setNames(cols[g1], cGenes)
#
#     #Set pseudogenes or ORF,... to light grey
#     ind <- which(gene.type$IMGT[[tsp]][names(TCRgene2aes[[tsp]][[s]]$color1)] %in% c("P", "ORF", "(ORF)"))
#     TCRgene2aes[[tsp]][[s]]$color1[ind] <- "grey95"
#
#     TCRgene2aes[[tsp]][[s]]$shape1 <- setNames(shapes[g2], cGenes)
#     TCRgene2aes[[tsp]][[s]]$outerColor1 <- setNames(cols_out[g2], cGenes)
#   }
# }
# for(t in c("color1", "shape1", "outerColor1", "color2")){
#   TCRgene2aes$HomoSapiens$TRBV[[t]]["TRBV12-3/12-4"] <- TCRgene2aes$HomoSapiens$TRBV[[t]]["TRBV12-3"]
#   TCRgene2aes$HomoSapiens$TRBV[[t]]["TRBV6-2/6-3"] <- TCRgene2aes$HomoSapiens$TRBV[[t]]["TRBV6-2"]
# }
#
# # If we want to make figures showing the color/shape from each gene.
# if (TRUE){
#   for (cScheme in 1:2){
#     for (tsp in species.list){
#       gg <- list()
#       for (s in segment.list){
#         cGenes <- names(TCRgene2aes[[tsp]][[s]][[paste0("color", cScheme)]])
#         cGenes <- setdiff(cGenes, "Other")
#         # Won't show the color used for "Other" in this plot as doesn't
#         # correspond to a gene.
#         nGenes <- length(cGenes)
#         cTab <- data.frame(gene = cGenes,
#           x = rep(1:3, each=ceiling(nGenes/3))[1:nGenes],
#           y = -rep(1:ceiling(nGenes/3), length.out=nGenes))
#         colorScale <- TCRgene2aes[[tsp]][[s]][[paste0("color", cScheme)]]
#         if (paste0("shape", cScheme) %in% names(TCRgene2aes[[tsp]][[s]])){
#           shapeScale <- TCRgene2aes[[tsp]][[s]][[paste0("shape", cScheme)]]
#         } else {
#           shapeScale <- rep(21, nGenes)
#         }
#         if (paste0("outerColor", cScheme) %in% names(TCRgene2aes[[tsp]][[s]])){
#           outerColorScale <- TCRgene2aes[[tsp]][[s]][[paste0("outerColor", cScheme)]]
#         } else {
#           outerColorScale <- rep("gray20", nGenes)
#         }
#         gg[[s]] <- ggplot(cTab, aes(x, y)) +
#           geom_point(aes(shape = gene, fill=gene, color=gene), stroke=1.5,
#             size=ifelse(nGenes <= 60, 6, 5)) +
#           geom_text(aes(label = gene), hjust = 0, nudge_x = 0.15) +
#           scale_fill_manual(values=colorScale) +
#           scale_shape_manual(values=shapeScale) +
#           scale_color_manual(values=outerColorScale) +
#           scale_x_continuous(limits=c(1,3.5)) +
#           theme_void() +
#           theme(legend.position = "none") +
#           ggtitle(s) +
#           theme(plot.title = element_text(size = 16, hjust=0.5, face="bold"),
#             plot.background=element_rect(color="gray30"))
#       }
#       fig <- ggpubr::ggarrange(plotlist=gg, nrow=2, ncol=2)
#       fig <- ggpubr::annotate_figure(fig,
#         top = ggpubr::text_grob(paste0(tsp, " - colorScheme ", cScheme),
#           color = "black", face = "bold", size = 16))
#       dir <- paste("figures/", sep="")
#       if (!dir.exists(dir)){
#         dir.create(dir)
#       }
#       filename=paste0(dir, "TCR_genes_tables_scheme", cScheme, "_", tsp, ".pdf")
#       ggsave(fig, filename=filename, device="pdf",
#         width = 15, height = 20)
#     }
#   }
# }
#
# # Defining other parameters -----------------------------------------------
# # This initiates a few values and load the mapping of CDR1/2 sequences from
# # the gene*allele names.
#
# Lmin <- 7
# Lmax <- 22
#
# aa <- "ACDEFGHIKLMNPQRSTVWY"
# aa.list <- unlist(strsplit(aa, split=""))
# N.aa <- length(aa.list)
# #map.aa <- 1:20
# #names(map.aa) <- aa.list
#
# gap <- c("-", "X", ".", "*","_")  # "x" is kept for 'real' gaps in CDR1/2, while all the other symbols are for missing data.
# th <- theme(plot.title = element_text(size = 8, hjust=0.5),
#   axis.title=element_text(size=8))
# yl <- ylim(0,2)
#
# chain.list <- c("TRA", "TRB")
# chain.small <- c("a","b")
# names(chain.small) <- chain.list
#
# species.list <- c("HomoSapiens", "MusMusculus")
#
# ### Load the aligned CDR1, CDR2 and CDR3 sequences
#
# cdr123 <- list()
# Jseq <- list()
#
# for(species in species.list){
#
#   tcdr123 <- list()
#   tJseq <- list()
#
#   for(chain in chain.list){
#     tcdr123[[chain]] <- read.csv(file=paste("data_raw/CDR123/",species,"/",chain,"V_allele.csv", sep=""), row.names = 1)
#     tJseq[[chain]] <- read.csv(file=paste("data_raw/CDR123/",species,"/",chain,"J_allele.csv", sep=""), row.names = 1)
#
#     #Take the data without allele information
#     tcdr123[[chain]] <- rbind(tcdr123[[chain]],read.csv(file=paste("data_raw/CDR123/",species,"/",chain,"V.csv", sep=""), row.names = 1))
#     tJseq[[chain]] <- rbind(tJseq[[chain]],read.csv(file=paste("data_raw/CDR123/",species,"/",chain,"J.csv", sep=""), row.names = 1))
#
#     #Change the gap ("-") int "x"
#     tcdr123[[chain]][] <- data.frame(lapply(tcdr123[[chain]], function(x){gsub(pattern="-", replacement = "g", x)}))
#     tJseq[[chain]][] <- data.frame(lapply(tJseq[[chain]], function(x){gsub(pattern="-", replacement = "g", x)}))
#
#
#   }
#
#   cdr123[[species]] <- list(tcdr123[["TRA"]],tcdr123[["TRB"]])
#   names(cdr123[[species]]) <- c("TRA","TRB")
#   Jseq[[species]] <- list(tJseq[["TRA"]],tJseq[["TRB"]])
#   names(Jseq[[species]]) <- c("TRA","TRB")
#
# }
#
# #Add CDR3 for TRBV6-2/6-3 and TRBV12-3/12-4
# rn <- rownames(cdr123$HomoSapiens$TRB)
# v6 <- cdr123$HomoSapiens$TRB["TRBV6-2",]
# names(v6) <- colnames(cdr123$HomoSapiens$TRB)
# v12 <- c("","",cdr123$HomoSapiens$TRB["TRBV12-3",c("CDR3", "default_allele","gene_type")],"", "")
# names(v12) <- colnames(cdr123$HomoSapiens$TRB)
# cdr123$HomoSapiens$TRB <- rbind(cdr123$HomoSapiens$TRB,v6,v6,v12,v12)
# rownames(cdr123$HomoSapiens$TRB) <- c(rn,"TRBV6-2/6-3", "TRBV6-2/6-3*01","TRBV12-3/12-4", "TRBV12-3/12-4*01")
#
#
# #Build reference sequences for the beginning and end of CDR3, allowing mapping to multiple alleles in most cases
#
# ref.cdr3.first <- list()
# ref.cdr3.last <- list()
# for(species in species.list){
#
#   ref.cdr3.first[[species]] <- list()
#   ref.cdr3.last[[species]] <- list()
#
#   for(chain in c("TRA", "TRB")){
#
#     V.end <- cdr123[[species]][[chain]][,"CDR3"]
#     V.list <- rownames(cdr123[[species]][[chain]])
#     V.list2 <- sapply(V.list, function(x){unlist(strsplit(x,split="*", fixed=T))[1]})
#     ref.cdr3.first[[species]][[chain]] <- lapply(V.list, function(x){
#       gn <- unlist(strsplit(x,split="*", fixed=T))[1]
#       ind <- which(V.list2==gn & V.end != "")
#       ref <- unique(V.end[ind])
#       #Remove cases which are subset of a longer one
#       lg.max <- max(nchar(ref))
#       ref.max <- ref[nchar(ref)==lg.max]
#       ref.final <- ref.max
#       for(r in ref){
#         if(nchar(r) != lg.max){
#           if(!r %in% substr(ref.max,1,nchar(r))) {
#             ref.final <- c(ref.final,r)
#           }
#         }
#       }
#       return(ref.final)
#     })
#
#     names(ref.cdr3.first[[species]][[chain]]) <- V.list
#
#
#     J.start <- Jseq[[species]][[chain]][,"CDR3"]
#     J.list <- rownames(Jseq[[species]][[chain]])
#     J.list2 <- sapply(J.list, function(x){unlist(strsplit(x,split="*", fixed=T))[1]})
#     ref.cdr3.last[[species]][[chain]] <- lapply(J.list, function(x){
#       gn <- unlist(strsplit(x,split="*", fixed=T))[1]
#       ind <- which(J.list2==gn & J.start != "")
#       ref <- unique(J.start[ind])
#       return(ref)
#     })
#     names(ref.cdr3.last[[species]][[chain]]) <- J.list
#   }
# }
#
#
# # Print the cases with multiple references
# print.ref <- F
# if(print.ref){
#   for(species in species.list){
#     for(chain in c("TRA", "TRB")){
#       print(c(species,chain))
#
#       lg <- sapply(ref.cdr3.first[[species]][[chain]],length)
#       for(i in which(lg!=1 & !grepl("*", names(lg), fixed=T))){
#         print(ref.cdr3.first[[species]][[chain]][i])
#       }
#       lg <- sapply(ref.cdr3.last[[species]][[chain]],length)
#       for(i in which(lg!=1 & !grepl("*", names(lg), fixed=T))){
#         print(ref.cdr3.last[[species]][[chain]][i])
#       }
#     }
#   }
# }

### Load the Table for flagging sequences

T_V <- list()
T_J <- list()

for(species in c("HomoSapiens", "MusMusculus")){

  temp_T_V <- list()
  temp_T_J <- list()

  for(chain in c("TRA", "TRB")){
    temp_T_V[[chain]] <- read.csv(file=paste("data_raw/Tables_T/",species,"/df_",chain,"V_substrings_to_check.csv", sep=""), row.names = 1)
    temp_T_J[[chain]] <- read.csv(file=paste("data_raw/Tables_T/",species,"/df_",chain,"J_substrings_to_check.csv", sep=""), row.names = 1)

  }

  T_V[[species]] <- list(temp_T_V[["TRA"]],temp_T_V[["TRB"]])
  names(T_V[[species]]) <- c("TRA","TRB")
  T_J[[species]] <- list(temp_T_J[["TRA"]],temp_T_J[["TRB"]])
  names(T_J[[species]]) <- c("TRA","TRB")

}


#
# #Establish a mapping of alternative column names
#
# mapping.colnames <- list()
#
# m <- read.csv("data_raw/TidyVJ/mapping_colnames.csv")
# for(chain in c("A","B","AB")){
#   ind <- which(grepl(chain,m[,"chain"]))
#   mapping.colnames[[chain]] <- m[ind,"CorrectedName"]
#   names(mapping.colnames[[chain]]) <- m[ind,"Name"]
# }
#
# #Data for formats based on clone.id
#
# clone.format.col <- list()
# clone.format.col[["AIRR"]] <- c("v_call", "j_call", "junction_aa")
# clone.format.col[["10X"]] <- c("v_gene", "j_gene", "cdr3")
# clone.format.col[["MiXCR"]] <- c("allVHitsWithScore", "allJHitsWithScore", "aaSeqCDR3")
# clone.format.col[["Qiagen"]] <- c("V-region", "J-region", "CDR3 amino acid seq")
# clone.format.col[["Adaptive.v4"]] <- c("v_resolved", "j_resolved", "amino_acid")
# clone.format.col[["Adaptive"]] <- c("vGeneName", "jGeneName", "aminoAcid")
# clone.format.col[["VDJdb"]] <- c("V", "J", "CDR3")
#
# clone.id <- c("clone_id", "cell_id", "cloneId", "barcode", "complex.id")
#
# test_file <- read.csv("data_raw/test.csv")

# Saving all these variable for internal use within the package ------------------
# functions

# usethis::use_data(gene.allele.list, gene.list, allele.default,
#   merge.mouse.TRAV, map, cdr123, Jseq, ref.cdr3.first, ref.cdr3.last,
#   th, yl, aa, aa.list, N.aa, mapping.colnames, clone.format.col, clone.id,
#   chain.small, gap, Lmin, Lmax, species.list, TCRgene2aes, T_V, T_J, test_file,
#   overwrite=T, internal=T)
usethis::use_data( T_V, T_J,
                   overwrite=T, internal=T)
usethis::use_data( T_V, T_J,
                  overwrite=T, internal=F)
# usethis::use_data(gene.allele.list, gene.list, allele.default, clone.format.col, clone.id,
#   merge.mouse.TRAV, map, cdr123, Jseq, ref.cdr3.first, ref.cdr3.last, mapping.colnames, T_V, T_J, test_file,
#   overwrite=T, internal=F)


# Resaving also the baseline_xxx.rds file for simple loading within -------------
# the package. Instead of having these files saved as .rds files and then loaded
# when needed through all.baseline <- readRDS(...), I read these here will save
# them as rda variables as well so that we can then directly load them
# through "all.baseline <- MixTCRviz::baseline_xxx" (otherwise path to these data
# isn't necessarily known).
# If we don't want that these are available for every use outside of the package
# but only to be used from within the package, we could save them as internal
# variables (but it's slowing the installation of the package and its loading
# due to the size of these data). We could also put these in "inst/RData/"
# folder for example and then in the code, use
# "system.file("RData", package = "MixTCRviz")" to get path to this folder.

