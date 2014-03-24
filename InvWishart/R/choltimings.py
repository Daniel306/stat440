
import matplotlib.pyplot as plt
import pandas
import scipy

R = pandas.read_csv("randcov.times.Rdata")

C = 5e9 #scale factor, to make things legible
#C = 5e6 #a scale factor which, if you turn off d^3, demonstrates (I think?) that this curve cannot be d^2: it initially lies below d^2 and then goes above it.
#C /= 100
print(C)
# C = 5e9 seems to fit the data well
#C *= 222
plt.plot(R["X1"], R["X1"]**2/C, label="d^2")
plt.plot(R["X1"], R["X1"]**3/C, label="d^3")
#plt.plot(R["X1"], R["X1"]**4/C, label="d^4")
#plt.plot(R["X1"], scipy.exp(R["X1"])/C, label="e^d")
plt.scatter(R["X1"], R["X0"])
plt.legend()

plt.show()
