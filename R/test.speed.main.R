# test.speed.main.R
#
#  Until we learn the proper Rish way of writing tests.
#  This script will serve.
#  in Python, the contents of this file would be under the
#   "if __name__ == '__main__'" 
#  block.

source("test.speed.R")


main <- function() { #avoid polluting the namespace
    
    results = NIW.runtimes(c("extremelynaive", "naive", "snappy1", "snappy2", "snappy3"), max=1e3, rep=1)

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

}

main()