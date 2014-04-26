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
    results = NIW.runtimes(krNIWAlgorithms$algorithm, max=50000, freq=3, rep=1) # Thorough: max=1e6, rep=5)

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
    for(int.a in 1:length(krNIWAlgorithms$algorithm)) {
        a = krNIWAlgorithms[int.a,]$algorithm
        a.results = results[results$algorithm == a, ] #don't forget the comma!
        lines(a.results$n, a.results$time, col = as.character(krNIWAlgorithms$color[int.a]), lty= as.character(krNIWAlgorithms$linestyle[int.a]))
        points(a.results$n, a.results$time, col = "black")
    }
    message("                     ")

    qq = as.character(krNIWAlgorithms$color)
    #print(qq) #DEBUG
    yy = as.character(krNIWAlgorithms$linestyle)
    #print(yy) #DEBUG
    legend("topleft", as.character(krNIWAlgorithms$algorithm), col=qq, lty=yy)
    
    message("Summary table:")
    message("(slope is in samples/second, and ratio is the ratio of that; sd is the standard deviation of ratio)")
    table <- NIW.runtime.createTable(results,"naive") #extremelynaive")
    table = table[order(table$ratio), ] #sort by relative speed
    print(table)
}

main()

