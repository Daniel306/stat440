# test.speed.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.


message("Loading test harness")
source("test.speed.R")
source("test.constants.R")
message("                     ")

message("Loading extension module")
require("Rcpp")
Rcpp::sourceCpp("rNIW.cpp") #make sure the sourcing of this doesn't interfere with the runtime measurements
message("                     ")

main <- function() { #avoid polluting the namespace
    message("Collecting runtimes")
    results = NIW.runtimes(krNIWAlgorithms, max=50000, freq=3, rep=1) # Thorough: max=1e6, rep=5)

    message("Largest runtimes:") # DEBUG
    print(results[results$n == max(results$n),]) #DEBUG

    message("                     ")
    
    message("Plotting")
    # plot runtime plots yay
    pdf(file="runtimes.pdf")
    plot(c(), xlim=range(results$n), ylim=range(results$time),
        main="rNIW implementation runtimes",
        xlab="n",
        ylab="time (s)")

    # TODO: can this whole loop be replaced with a single lines() call?
    for(a in levels(results$algorithm)) {
        int.a = which(levels(results$algorithm) == a) #extract an integer id for the algorithm for use as a display style index
        a.results = results[results$algorithm == a, ] #don't forget the comma
        
        lines(a.results$n, a.results$time, col = krNIWAlgorithms.colors[int.a], lty= krNIWAlgorithms.linestyle[int.a])
    }
    message("                     ")
   
    legend("topleft", levels(results$algorithm), col=krNIWAlgorithms.colors, lty=krNIWAlgorithms.linestyle)
    
    message("Summary table:")
    message("(slope is in samples/second, and ratio is the ratio of that; sd is the standard deviation of ratio)")
    table <- NIW.runtime.createTable(results,"naive") #extremelynaive")
    table = table[order(table$ratio), ] #sort by relative speed
    print(table)
}

main()

