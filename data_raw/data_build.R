# Description ---------------------------------------------------------------
# In this file we create and export some variables used by the package
# The script should be run from MixTCRviz folder
library(MixTCRviz)
### Load the Table for flagging sequences

T_V <- list()
T_J <- list()

for(species in c("HomoSapiens", "MusMusculus")){

  temp_T_V <- list()
  temp_T_J <- list()

  for(chain in c("TRA", "TRB")){
    temp_T_V[[chain]] <- read.csv(file=paste("data_raw/Tables_T/",species,"/df_",chain,"V_substrings_to_check.csv", sep=""), row.names = 1)
    temp_T_J[[chain]] <- read.csv(file=paste("data_raw/Tables_T/",species,"/df_",chain,"J_substrings_to_check.csv", sep=""), row.names = 1)

    names(temp_T_V[[chain]]) <- sub("^len_", "L_", names(temp_T_V[[chain]]))
    names(temp_T_J[[chain]]) <- sub("^len_", "L_", names(temp_T_J[[chain]]))
  }

  T_V[[species]] <- list(temp_T_V[["TRA"]],temp_T_V[["TRB"]])
  names(T_V[[species]]) <- c("TRA","TRB")
  T_J[[species]] <- list(temp_T_J[["TRA"]],temp_T_J[["TRB"]])
  names(T_J[[species]]) <- c("TRA","TRB")

}


## Load species.list from MixTCRviz
species.list <- MixTCRviz:::species.list
correct.VJnames <- MixTCRviz:::correct.VJnames
verify.chain <- MixTCRviz:::verify.chain

usethis::use_data( species.list, correct.VJnames, verify.chain,
                   overwrite=T, internal=T)
usethis::use_data( T_V, T_J,
                   overwrite=T, internal=F)
