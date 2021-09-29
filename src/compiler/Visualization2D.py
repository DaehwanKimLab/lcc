import csv
import matplotlib.pyplot as plt
import numpy as np

def VisualizeData(Dataset, RequestedData):
    if RequestedData['ID_Molecules']:
        PlotData(Dataset, RequestedData['ID_Molecules'])
    if RequestedData['Name_Variables']:
        PlotData(Dataset, RequestedData['Name_Variables'])


def PlotData(Dataset, TargetNames):

    X = list()
    Y = list()

    with open(Dataset) as CsvFile:
        Reader = csv.reader(CsvFile, delimiter=',')

        Name_Row = next(Reader)
        ColumnsForTargetNames = []
        for TargetName in TargetNames:
            ColumnsForTargetNames.append(Name_Row.index(TargetName))
            # print("found a column for %s" % TargetName)

        NumRows = 0
        for Row in Reader:
            ValuesInColumnsForTargetNames = list()
            for Column in ColumnsForTargetNames:
                ValuesInColumnsForTargetNames.append(int(Row[Column]))
            Y.append(ValuesInColumnsForTargetNames)
            NumRows += 1

        # print("Rows = ", numRows)
        X = list(range(NumRows))
        # print(len(X) == len(Y))

    Y = np.array(Y).transpose()

    # ax = plt.axes(xlim=(0, len(X)), ylim=(0, Y.max() * 1.2))
    ax = plt.axes()

    ax.set_xbound(0, len(X))
    ax.set_ybound(0, Y.max() * 1.2)
    ax.set_ylabel('Count')
    ax.set_xlabel('Time (seconds) ')
    ax.set_title(TargetNames)

    for i, TargetName in enumerate(TargetNames):
        ax.plot(X, Y[i], label=TargetName)

    ax.legend(loc='upper left')

    plt.show()


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