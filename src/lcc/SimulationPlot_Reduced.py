import matplotlib.pyplot as plt

def MMequation(Vmax, S, km):
    return (Vmax * S) / (km + S)

# Concentration
A = 50
B = 0
Vmax = 25
km = 10

TimeResolution = 100

Data_A = list()
Data_B = list()
Data_Rate = list()

Data_A.append(A)
Data_B.append(B)
Data_Rate.append(0)

SimSteps = 5000

i = 0
while i < SimSteps:

    Rate = MMequation(Vmax, A, 13) / TimeResolution

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

ax2.plot(X, Y, color="blue")
ax2.set_title('Phase plane')
ax2.set_xlabel("[A]")
ax2.set_ylabel("[B]")
ax2.grid()

ax3.plot(X, 'r-', label="[A]")
ax3.plot(R, 'b-', label="Rate")
ax3.set_title('Rate vs [A]')
ax3.set_xlabel('[A]')
ax3.legend(loc='upper left')
ax3.grid()

plt.show()