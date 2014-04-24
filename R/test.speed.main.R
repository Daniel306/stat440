# test.speed.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.

source("test.speed.R")
source("test.constants.R")

main <- function() { #avoid polluting the namespace

  # TODO: pull the alg names out into a constant..    
    results = NIW.runtimes(krNIWAlgorithms, max=1e6, rep=3)

    # plot runtime plots yay
    plot(c(), xlim=range(results$n), ylim=range(results$time),
        main="rNIW implementation runtimes",
        xlab="n",
        ylab="time (s)")
    for(a in levels(results$algorithm)) {
        int.a = which(levels(results$algorithm) == a) #extract an integer id for the algorithm for use as a colour palette index
        a.results = results[results$algorithm == a, ] #don't forget the comma
        
        lines(a.results$n, a.results$time, col = int.a, lty="solid")
    }

    qq=1:length(levels(results$algorithm))
    #print(qq) #DEBUG
    yy = rep("solid", length(levels(results$algorithm)))
    #print(yy) #DEBUG
    legend("topleft", levels(results$algorithm), col=qq, lty=yy)

    table <- NIW.runtime.createTable(results,"extremelynaive")
    print(table)
}

main()



#example from one run
#algorithm     ratio        sd
#1     Rcpp2 8.7011074 2.5359325
#2   snappy2 0.7468129 0.1213526
#3     naive 1.0000000 0.1674202