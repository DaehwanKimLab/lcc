import os
import csv
import matplotlib.pyplot as plt
import numpy as np

def LoadRawData(Data_Dir):
    Datasets = dict()

    def Parse_tsv(FilePath, FileName):
        Fullpath = FilePath + '/' + FileName
        # print(fname)

        with open(Fullpath) as fp:
            Reader = csv.reader(fp, delimiter=',')
            Datasets[FileName] = list(Reader)

    for FileName in os.listdir(Data_Dir):
        if FileName.endswith('.tsv'):
            Parse_tsv(Data_Dir, FileName)

    return Datasets

def PlotData(Dataset):

    Legend = Dataset[0]
    Data = Dataset[1:]

    for i, Row in enumerate(Data):
        Data[i] = [int(NumStr) for NumStr in Row]

    Idx_SimStep = 1
    Idx_Vol = 1 + Idx_SimStep
    Idx_Enzyme = 8 + Idx_Vol
    Idx_Substrate = 21 + Idx_Enzyme

    Data_Transposed = np.array(Data).transpose()
    Data_SimStep = Data_Transposed[0:Idx_SimStep]
    Data_Vol = Data_Transposed[Idx_SimStep:Idx_Vol]
    Data_EnzCounts = Data_Transposed[Idx_Vol:Idx_Enzyme]
    Data_SubCounts = Data_Transposed[Idx_Enzyme:Idx_Substrate]

    Legend_SimStep = Legend[0:Idx_SimStep]
    Legend_Vol = Legend[Idx_SimStep:Idx_Vol]
    Legend_EnzCounts = Legend[Idx_Vol:Idx_Enzyme]
    Legend_SubCounts = Legend[Idx_Enzyme:Idx_Substrate]

# Potentially split point for another function

    X = Data_SimStep
    Y = Data_SubCounts

    assert X.shape[1] == Y.shape[1]

    ax = plt.axes(xlim=(0, len(X)), ylim=(0, Y.max() * 1.2))

    ax.set_ylabel('Substrate Count')
    ax.set_xlabel('Sim Step')
    ax.set_title('TCA cycle')

    for i in range(len(Legend_SubCounts)):
        ax.plot(X[0], Y[i], label=Legend_SubCounts[i])

    ax.legend(loc='upper left')

    plt.show()

Data_Dir = '../lcc'
Datasets = LoadRawData(Data_Dir)
for FileName, Dataset in Datasets.items():
    PlotData(Dataset)


# SavePath = 'lccsave/'
# FileName = 'DL_EcoliSimulation_20210928-1517_Transcription.csv'
# Dataset = SavePath + FileName
#
# RequestedData_Test = dict()
# RequestedData_Test['ID_Molecules'] = ['nusB', 'ftsZ']
# RequestedData_Test['Name_Variables'] = None
#
# VisualizeData(Dataset, RequestedData_Test)



'''

# Ronnin's Original Plot Animation Code
 
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

def rect(x, B):
    """
    create a rectangle function
    returns a numpy array that is 1 if |x| < w and 0 if |x| > w
    B is the rectangle width centered at 0
    x is the number of points in the array
    """

    B = int(B)
    x = int(x)

    high = np.ones(B)
    low1 = np.zeros(int(x / 2 - B / 2))
    x1 = np.append(low1, high)
    rect = np.append(x1, low1)

    if x > len(rect):
        rect = np.append(rect, 0)
    elif x < len(rect):
        rect = rect[:-1]

    return rect


maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)

x = []
y = []
targetname = 'nusB'
with open('DL_EcoliSimulation_20210928-1517_Transcription.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')

    name_row = next(reader)
    colNum = 0
    for name in name_row:
        if name == targetname:
            print("found GLB")
            break
        colNum += 1

    numRows = 0
    for row in reader:
        y.append(int(row[colNum]))
        numRows += 1

    print("rows = ", numRows)
    x = list(range(numRows))
    print(len(x) == len(y))
    # print(x)
    # print(y)




# %matplotlib notebook


fig = plt.figure()
ax = plt.axes(xlim=(0, 600), ylim=(0, 10))

ax.set_ylabel(targetname + ' Count')
ax.set_xlabel('Time (seconds) ')
ax.set_title(targetname + ' Count vs. Time')
lines = plt.plot([], 'b-', markersize=2)
line = lines[0]

xanim = []
yanim = []

# set data and just scroll through xlims..?
line.set_data((x, y))
offset = 850


def animate(i):
    index = i + offset
    if index < len(x):
        xanim.append(x[index])
        yanim.append(y[index])
        ax.set_xlim(min(xanim), index + 10)
        ax.set_ylim(min(yanim) - 10, max(yanim) + 10)

        line.set_data((xanim, yanim))


ani = FuncAnimation(fig, animate, frames=2000, interval=500)
plt.show()

f = r"test.gif"
writergif = animation.PillowWriter(fps=30)
ani.save(f, writer=writergif)
'''