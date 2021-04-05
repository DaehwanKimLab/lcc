import matplotlib.pyplot as plt

def VisualizeData(Title, X, Y):
    fig, ax = plt.subplots()
    ax.plot(X, Y)

    ax.set(xlabel='time (sim steps)', ylabel='con', title=Title)
    ax.grid()

    fig.savefig(Title)

    return(plt.show())
