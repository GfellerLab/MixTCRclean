library(MixTCRviz)

### Load the Table for flagging sequences
T_V <- list()
T_J <- list()

thresholds <- c("0_95", "0_98", "0_99", "0_995")

for(thr in thresholds){
  temp_T_V_species <- list()
  temp_T_J_species <- list()

  for(species in c("HomoSapiens", "MusMusculus")){
    temp_T_V <- list()
    temp_T_J <- list()

    for(chain in c("TRA", "TRB")){
      temp_T_V[[chain]] <- read.csv(file=paste0("data_raw/Tables_T/", species, "/df_", chain, "V_substrings_to_match_symm_window_width_3_thr_", thr, ".csv"), row.names = 1)
      temp_T_J[[chain]] <- read.csv(file=paste0("data_raw/Tables_T/", species, "/df_", chain, "J_substrings_to_match_symm_window_width_3_thr_", thr, ".csv"), row.names = 1)
    }

    temp_T_V_species[[species]] <- list(TRA = temp_T_V[["TRA"]], TRB = temp_T_V[["TRB"]])
    temp_T_J_species[[species]] <- list(TRA = temp_T_J[["TRA"]], TRB = temp_T_J[["TRB"]])
  }

  T_V[[thr]] <- temp_T_V_species
  T_J[[thr]] <- temp_T_J_species
}

## Load species.list from MixTCRviz
species.list <- MixTCRviz:::species.list
correct.VJnames <- MixTCRviz:::correct.VJnames
verify.chain <- MixTCRviz:::verify.chain
usethis::use_data(species.list, correct.VJnames, verify.chain,
                  overwrite = TRUE, internal = TRUE)
usethis::use_data(T_V, T_J,
                  overwrite = TRUE, internal = FALSE)
