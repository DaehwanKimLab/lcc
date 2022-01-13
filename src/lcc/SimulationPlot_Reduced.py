import matplotlib.pyplot as plt

def MMequation(Vmax, S, KM):
    return (Vmax * S) / (KM + S)

# Concentration
A = 50
B = 0
Vmax = 50
KM = 15

Data_A = list()
Data_B = list()
Data_Rate = list()

Data_A.append(A)
Data_B.append(B)
Data_Rate.append(0)

SimSteps = 5000
TimeResolution = 1000

i = 0
while i < SimSteps:

    Rate = MMequation(Vmax, A, KM) / TimeResolution

    A = A - Rate
    B = B + Rate

    Data_A.append(A)
    Data_B.append(B)
    Data_Rate.append(Rate)

    i += 1

X = Data_A
Y = Data_B
R = Data_Rate

fig = plt.figure()
fig.subplots_adjust(wspace=0.2, hspace=0.3)
# Dynamics
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

ax1.plot(X, 'r-', label="[A]")
ax1.plot(Y, 'b-', label="[B]")
ax1.set_title('Dynamics')
ax1.set_xlabel('SimStep')
ax1.legend(loc='upper left')
ax1.grid()

# ax = plt.axes(xlim=(0, X.max()), ylim=(0, Y.max()))
ax2.plot(X, 'r-', label="[A]")
ax2.plot(R, 'b-', label="Rate")
ax2.set_title('Rate vs [A]')
ax2.set_xlabel('[A]')
ax2.legend(loc='upper left')
ax2.grid()

# Phase plane
ax3.plot(X, Y, color="blue")
ax3.set_title('Phase plane')
ax3.set_xlabel("[A]")
ax3.set_ylabel("[B]")
ax3.grid()

plt.show()