# test.speed.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.

source("test.speed.R")


main <- function() { #avoid polluting the namespace
    
    NIW.runtimes(c("snappy1", "snappy2"), max=1e3)


}

main()