#Run it from the MixTCRclean/ folder

#Do this if you do not want to install the package
devtools::load_all(".")

#Do this if you have installed the package
#library(MixTCRclean)


MixTCRclean(input=MixTCRclean::test_file, filename.output = "test",
            output.path="test/out/", verbose = 2, check.cdr3.mode = 2,
            correct.gene.names = T)

# Or if you downloaded the github repository:
# MixTCRclean(input="test/test.csv",
#             output.path="test/out/", verbose = 2, check.cdr3.mode = 2,
#             correct.gene.names = T)

help("MixTCRclean")
new.data <- F

if(!new.data){
  m <- read.csv("test/out/test_processed.csv")
  m.comp <- read.csv("test/out_compare/test_processed.csv")

  comp <- identical(m, m.comp) # should be TRUE
  if(comp){
    print("No problem detected")
  } else {
    stop("There were some issues...")
  }
} else {
  system("cp test/out/test_processed.csv test/out_compare/test_processed.csv")
}
