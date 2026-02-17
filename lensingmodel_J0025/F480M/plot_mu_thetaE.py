import os, sys
import numpy as np
import matplotlib.pyplot as plt

sigma, mu, thetaE = np.loadtxt('muthetaE.txt')


# Using this plot to find the maximum Einstein radius and the corresponding magnification.
plt.plot(thetaE, mu)
plt.show()



