#Run it from the MixTCRclean/ folder

#Do this if you do not want to install the package
devtools::load_all(".")

#Do this if you have installed the package
#library(MixTCRclean)



MixTCRclean(input1="test/test.csv",
            output.path="test/out/", verbose = 2, check.cdr3.mode = 2,
            correct.gene.names = T)

new.data <- F

if(!new.data){
  m <- read.csv("test/out/processed_data/test.csv")
  m.comp <- read.csv("test/out_compare/test.csv")

  comp <- identical(m, m.comp) # should be TRUE
  if(comp){
    print("No problem detected")
  } else {
    stop("There were some issues...")
  }
} else {
  system("cp test/out/processed_data/test.csv test/out_compare/test.csv")
}
