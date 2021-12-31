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
            Reader = csv.reader(fp, delimiter='\t')
            Datasets[FileName] = list(Reader)

    for FileName in os.listdir(Data_Dir):
        if FileName.endswith('.tsv'):
            Parse_tsv(Data_Dir, FileName)

    return Datasets

def PlotData(Dataset, Type):
    # Type 1 --> plot over sim step
    # Type 2 --> phase plane
    # Type 3 --> both

    Legend = Dataset[0]
    Data = Dataset[1:]

    for i, Row in enumerate(Data):
        Data[i] = [float(NumStr) for NumStr in Row]

    Data_Transposed = np.array(Data).transpose()
    Data_SimStep = Data_Transposed[0:1]

    # Default setting
    Title = ''

    plot_TCACycle = False
    plot_PolymeraseReactions = False
    plot_Ecoli_NoTCA = False
    plot_Ecoli_TCA = False
    plot_ShowAll = False

    # plot_TCACycle = True
    # plot_PolymeraseReactions = True
    # plot_Ecoli_NoTCA = True
    plot_Ecoli_TCA = True
    # plot_ShowAll = True

    # TCA Cycle
    if plot_TCACycle:
        Molecules = ['SimStep', 'Vol', 'GltA', 'oxaloacetate', 'acetyl-CoA', 'H2O', 'citrate', 'CoA', 'AcnA', 'isocitrate', 'Icd', 'NAD+', 'keto-glutarate', 'NADH', 'CO2', 'SucA', 'succinyl-CoA', 'SucD', 'ADP', 'Pi', 'succinate', 'ATP', 'Sdh', 'FAD', 'fumarate', 'FADH2', 'FumA', 'malate', 'Mdh', 'H+']

        Idx_SimStep = 1
        Idx_Vol = 1 + Idx_SimStep
        Idx_Enzyme = 8 + Idx_Vol

        Data_SimStep = Data_Transposed[0:Idx_SimStep]
        Data_Vol = Data_Transposed[Idx_SimStep:Idx_Vol]
        Data_EnzCounts = Data_Transposed[Idx_Vol:Idx_Enzyme]
        Data_SubCounts = Data_Transposed[Idx_Enzyme:]

        Legend_SimStep = Legend[0:Idx_SimStep]
        Legend_Vol = Legend[Idx_SimStep:Idx_Vol]
        Legend_EnzCounts = Legend[Idx_Vol:Idx_Enzyme]
        Legend_SubCounts = Legend[Idx_Enzyme:]

        Title = 'TCA cycle'
        Data = Data_SubCounts
        Legend = Legend_SubCounts

    # Polymerase Reactions
    elif plot_PolymeraseReactions:
        Idx_SimStep = 1
        Idx_Vol = 1 + Idx_SimStep
        Idx_Polymerase = 1 + Idx_Vol
        Idx_ReactionSubstrate = 1 + Idx_Polymerase

        Data_SimStep = Data_Transposed[0:Idx_SimStep]
        Data_Vol = Data_Transposed[Idx_SimStep:Idx_Vol]
        Data_PolymeraseCounts = Data_Transposed[Idx_Vol:Idx_Polymerase]
        Data_ReactionSubstrateCounts = Data_Transposed[Idx_Polymerase:Idx_ReactionSubstrate]
        Data_BuildingBlockCounts = Data_Transposed[Idx_ReactionSubstrate:]
        Data_AllSmallMoleculeCounts = Data_Transposed[Idx_Polymerase:]

        Legend_SimStep = Legend[0:Idx_SimStep]
        Legend_Vol = Legend[Idx_SimStep:Idx_Vol]
        Legend_PolymeraseCounts = Legend[Idx_Vol:Idx_Polymerase]
        Legend_ReactionSubstrateCounts = Legend[Idx_Polymerase:Idx_ReactionSubstrate]
        Legend_BuildingBlockCounts = Legend[Idx_ReactionSubstrate:]
        Legend_AllSmallMoleculeCounts = Legend[Idx_Polymerase:]

        Title = 'PolymeraseReaction'
        Data = Data_BuildingBlockCounts
        Legend = Legend_BuildingBlockCounts

    elif plot_Ecoli_NoTCA:
    # Temporary capping for plotting speed
        Title = 'Ecoli_noTCA'
        N_Species = 500
        GeneBegins = 35 + 2
        mRNABegins = 4610 + 2
        ProteinBegins = 9180 + 2
        # Target = GeneBegins
        # Target = mRNABegins
        Target = ProteinBegins

        Data = Data_Transposed[Target:Target+N_Species]
        Legend = Legend[Target:Target+N_Species]

    elif plot_Ecoli_TCA:

        # Reference for molecules before Ch1, Genes, RNAs, Proteins
        Molecules = ['SimStep', 'Vol', 'GltA', 'oxaloacetate', 'acetyl-CoA', 'H2O', 'citrate', 'CoA', 'AcnA', 'isocitrate',
                     'Icd', 'NAD+', 'keto-glutarate', 'NADH', 'CO2', 'SucA', 'succinyl-CoA', 'SucD', 'ADP', 'Pi',
                     'succinate', 'ATP', 'Sdh', 'FAD', 'fumarate', 'FADH2', 'FumA', 'malate', 'Mdh', 'H+', 'pol1', 'ppi',
                     'dATP', 'dCTP', 'dGTP', 'dUTP', 'rnap', 'CTP', 'GTP', 'UTP', 'r1', 'H_{2}O', 'ALA', 'ARG', 'ASN',
                     'ASP', 'CYS', 'GLT', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR',
                     'TRP', 'TYR', 'SEL', 'VAL']

        # Temporary capping for plotting speed
        Title = 'Ecoli_TCA'
        SmallMoleculesBegins = 0 + 2
        GeneBegins = 62 + 2
        mRNABegins = 4635 + 2
        ProteinBegins = 9176 + 2
        N_Species = 500

        # Target = GeneBegins
        # Target = mRNABegins
        # Target = ProteinBegins
        Target = SmallMoleculesBegins

        if Target == SmallMoleculesBegins:
            N_Species = GeneBegins - SmallMoleculesBegins - 1 # -1 For chromosome

        Data = Data_Transposed[Target:Target+N_Species]
        Legend = Legend[Target:Target+N_Species]

    # Show all
    else:
        Cap = 50
        Title = 'Show All'
        Data = Data_Transposed[2:Cap]
        Legend = Legend[2:Cap]


    # Plotting options
    if Type == 1:
        Plot_OverSimStep(Title, Data_SimStep[0], Data, Legend)

    elif Type == 2:
        Plot_PhasePlane(Title, Data_SimStep[0], Data, Legend)

    elif Type == 3:
        Plot_OverSimStep(Title, Data_SimStep[0], Data, Legend)
        Plot_PhasePlane(Title, Data_SimStep[0], Data, Legend)

# Potentially split point for another function

def Plot_OverSimStep(Title, SimStep, Data, Legend):
    X = SimStep
    Y = Data

    assert len(X) == Y.shape[-1]

    ax = plt.axes(xlim=(0, len(X)), ylim=(0, Y.max() * 1.2))

    ax.set_ylabel('Count')
    ax.set_xlabel('Sim Step')
    ax.set_title(Title + " over SimStep")

    for i in range(len(Legend)):
        ax.plot(X, Y[i], label=Legend[i])

    ax.legend(loc='upper left')

    plt.show()

def Plot_PhasePlane(Title, SimStep, Data, Legend):
    SinglePair = False
    AllPairs = False
    # SinglePair = True
    AllPairs = True

    # Pick two indices for data
    Molecules = ['GltA', 'oxaloacetate', 'acetyl-CoA', 'H2O', 'citrate', 'CoA', 'AcnA', 'isocitrate', 'Icd', 'NAD+', 'keto-glutarate', 'NADH', 'CO2', 'SucA', 'succinyl-CoA', 'SucD', 'ADP', 'Pi', 'succinate', 'ATP', 'Sdh', 'FAD', 'fumarate', 'FADH2', 'FumA', 'malate', 'Mdh', 'H+']

    SelectedMolecules = ['oxaloacetate', 'acetyl-CoA', 'citrate', 'isocitrate', 'keto-glutarate', 'succinyl-CoA', 'succinate', 'fumarate', 'malate']
    # Subset for convenient visualization, assuming acetyl-CoA level is constant
    SelectedMolecules = ['oxaloacetate', 'citrate', 'isocitrate', 'keto-glutarate', 'succinyl-CoA', 'succinate', 'fumarate', 'malate']

    Idx_SelectedMolecules = list()
    Dict_SelectedMoleculesIdx = dict()
    for SelectedMolecule in SelectedMolecules:
        i = 0
        for Molecule in Molecules:
            if Molecule == SelectedMolecule:
                Idx_SelectedMolecules.append(i)
                Dict_SelectedMoleculesIdx[SelectedMolecule] = i
                break
            i += 1

    if SinglePair:

        fig = plt.figure()

        i = Idx_SelectedMolecules[5]
        j = Idx_SelectedMolecules[4]

        X = Data[i]
        Y = Data[j]

        fig = plt.figure()
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)

        ax1.plot(X, 'r-', label=Molecules[i])
        ax1.plot(Y, 'b-', label=Molecules[j])
        ax1.set_title(Title + ' Dynamics')
        ax1.set_xlabel('SimStep')
        # ax2.legend(loc='upper left')
        ax1.grid()

        # ax = plt.axes(xlim=(0, X.max()), ylim=(0, Y.max()))

        ax2.plot(X, Y, color="blue")
        ax2.set_title(Title + ' Phase plane')
        ax2.set_xlabel(Molecules[i])
        ax2.set_ylabel(Molecules[j])
        ax2.grid()

        plt.show()

    # all pairs
    elif AllPairs:
        fig = plt.figure()
        fig.subplots_adjust(wspace=0.2, hspace=0.3)
        # PossiblePairs_Legend = [(a, b) for idx_a, a in zip(Idx_SelectedMolecules, SelectedMolecules) for idx_b, b in zip(Idx_SelectedMolecules, SelectedMolecules) if idx_a < idx_b]
        PossiblePairs_Idx = [(idx_a, idx_b) for i, [name_a, idx_a] in enumerate(Dict_SelectedMoleculesIdx.items()) for j, [name_b, idx_b] in enumerate(Dict_SelectedMoleculesIdx.items()) if (i + 1 == j)]
        PossiblePairs_Names = [(name_a, name_b) for i, [name_a, idx_a] in enumerate(Dict_SelectedMoleculesIdx.items()) for j, [name_b, idx_b] in enumerate(Dict_SelectedMoleculesIdx.items()) if (i + 1 == j)]

        N_PossiblePairs = len(PossiblePairs_Idx)
        n = 1
        for i, j in PossiblePairs_Idx:
            X = Data[i]
            Y = Data[j]

            ax1 = fig.add_subplot(2, int(N_PossiblePairs) + 1, 2 * n - 1)
            ax2 = fig.add_subplot(2, int(N_PossiblePairs) + 1, 2 * n)

            ax1.plot(X, 'r-', label=Molecules[i])
            ax1.plot(Y, 'b-', label=Molecules[j])
            ax1.set_title(Title + ' Dynamics')
            ax1.set_xlabel('SimStep')
            ax1.set_ylabel('Count')
            ax1.legend(loc='upper left')
            ax1.grid()

            # ax = plt.axes(xlim=(0, X.max()), ylim=(0, Y.max()))

            ax2.plot(X, Y, color="blue")
            ax2.set_title(Title + ' Phase plane')
            ax2.set_xlabel(Molecules[i])
            ax2.set_ylabel(Molecules[j])
            ax2.grid()

            n += 1

        plt.show()




Data_Dir = '.'
Datasets = LoadRawData(Data_Dir)
for FileName, Dataset in Datasets.items():
    PlotData(Dataset, 2)


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
