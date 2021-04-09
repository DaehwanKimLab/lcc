import matplotlib.pyplot as plt

# Load npy files for data plotting

# Load txt files for legend info parsing


def VisualizeData(X, Y, XLabel, YLabel, Legends, Title): # Make Labels, Legends, and Title optional
    fig, ax = plt.subplots()
    Lines = ax.plot(X, Y)
    ax.legends(Lines, Legends)
    ax.set(xlabel=XLabel, ylabel=YLabel, title=Title)
    ax.grid()

    fig.savefig(Title)

    return(plt.show())
