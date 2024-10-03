import matplotlib.pyplot as plt
import numpy as np

sig1 = [1 for i in range(100)]
sig2 = [0 for i in range(60)]

for i in range(10):
    sig1[i] = 0

a = 1/153.0899  # measured values

for i in range(60):
    sig2[i] = a*np.exp(-a*i*20)

conv = np.convolve(sig1, sig2)
convCom = [0 for i in range(159)]

for i in range(139):
    convCom[i+9] = max(conv)*(1-np.exp(-a*20*i))
#conv /= len(sig2) # Normalize

plt.figure()
plt.plot(sig2)
plt.figure()
plt.plot(conv)
plt.plot(convCom)

print(len(conv))
print(len(convCom))
print(conv-convCom)
