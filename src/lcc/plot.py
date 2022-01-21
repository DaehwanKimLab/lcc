import os, sys
from argparse import ArgumentParser
import csv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import numpy as np

SaveFilename = ''

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

def PlotData(Dataset):

    Legend = Dataset[0]
    Data = Dataset[1:]

    for i, Row in enumerate(Data):
        Data[i] = [float(NumStr) for NumStr in Row]

    Data_Transposed = np.array(Data).transpose()
    Data_Time = Data_Transposed[0:1]
    Data_Vol = Data_Transposed[1:2]

    ### Default setting ###
    Title = ''

    plot_TCACycle = False
    plot_TCACycle_Reduced = False
    plot_PolymeraseReactions = False
    plot_Ecoli_NoTCA = False
    plot_Ecoli_TCA = False
    plot_ShowAll_DynamicsOnly = False
    plot_ShowAll_DynamicsAndPhasePlain = False

    CheckRateExpected = False

    ### Cycle option ### 

    # plot_TCACycle = True
    # plot_TCACycle_Reduced = True
    plot_PolymeraseReactions = True
    # plot_Ecoli_NoTCA = True
    # plot_Ecoli_TCA = True
    # plot_ShowAll_DynamicsOnly = True
    # plot_ShowAll_DynamicsAndPhasePlain = True

    ### Debugging option ###
    CheckRateExpected = True


    # Detect TCA legend for appropriate display
    if "oxaloacetate" in Legend:
        plot_ShowAll_DynamicsOnly = False
        plot_Ecoli_TCA = True

    # Show all
    if plot_ShowAll_DynamicsOnly:
        Cap = 50
        Title = ''
        Data = Data_Transposed[2:Cap] / Data_Vol
        Legend = Legend[2:Cap]

        Plot_Dynamics(Title, Data_Time[0], Data, Legend)
        return 0

    if plot_ShowAll_DynamicsAndPhasePlain:
        Cap = 50
        Title = ''
        Data = Data_Transposed[2:Cap] / Data_Vol
        Legend = Legend[2:Cap]

        if 'RL' in Legend:
            Plot_Dynamics(Title, Data_Time[0], Data, Legend)

        else:
            Idx_Pairs = None
            Names_Pairs = None

            if 'oxaloacetate' in Legend:
                Idx_Pairs, Names_Pairs = Indexing(Legend, 'TCA')
            else:
                Idx_Pairs, Names_Pairs = Indexing(Legend)

            Plot_DynamicsAndPhasePlaneAndRate(Title, Data_Time[0], Data, Legend, Idx_Pairs, Names_Pairs)

        return 0

    # Polymerase Reactions
    if plot_PolymeraseReactions:
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
        
        Plot_Dynamics(Title, Data_SimStep[0], Data, Legend)
        return 0

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

    elif plot_TCACycle_Reduced:
        Molecules = ['SimStep', 'Vol', 'GltA', 'oxaloacetate', 'acetyl-CoA', 'H2O', 'citrate', 'CoA']

        Idx_SimStep = 1
        Idx_Vol = 1 + Idx_SimStep
        Idx_Enzyme = 1 + Idx_Vol

        Data_SimStep = Data_Transposed[0:Idx_SimStep]
        Data_Vol = Data_Transposed[Idx_SimStep:Idx_Vol]
        Data_EnzCounts = Data_Transposed[Idx_Vol:Idx_Enzyme]
        Data_SubCounts = Data_Transposed[Idx_Enzyme:]

        Legend_SimStep = Legend[0:Idx_SimStep]
        Legend_Vol = Legend[Idx_SimStep:Idx_Vol]
        Legend_EnzCounts = Legend[Idx_Vol:Idx_Enzyme]
        Legend_SubCounts = Legend[Idx_Enzyme:]

        Title = 'TCA cycle_Reduced'
        Data = Data_SubCounts
        Legend = Legend_SubCounts

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


    # All TCA Cases
    Idx_Pairs, Names_Pairs = Indexing(Legend, Query="TCA")

    if CheckRateExpected:
        RateCheck(Data, Legend, Idx_Pairs, Names_Pairs)

    Plot_DynamicsAndPhasePlaneAndRate(Title, Data_Time[0], Data, Legend, Idx_Pairs, Names_Pairs)


def RateCheck(Data, Legend, Idx_Pairs, Names_Pairs):
    return 0
    for a, b in Idx_Pairs:
        Conc_A = Data[a]
        Conc_B = Data[b]
        assert len(Conc_A) == len(Conc_B)

        Conc_Enz = 0
        kcat = 0
        KM = 0
        if ('oxaloacetate' in Names_Pairs) and ('citrate' in Names_Pairs):
            Idx_Enz = Legend.index('GltA')  # TODO:make a dictionary
            Conc_Enz = Data[Idx_Enz]
            # kcat =
            # KM =
        else:
            print('No matching Conc_Enz is not determined')


        for i in range(len(Conc_A)):
            Conc_A_Before = Conc_A[i]
            Conc_A_After = Conc_A[i + 1]

            Conc_B_Before = Conc_B[i]
            Conc_B_After = Conc_B[i + 1]

            # Rate = MichaelisMentenEqn(Conc_Enzyme, Conc_Substrate, kcat, KM)
            # assert Conc_oxaloacetate - Conc_oxaloacetate_Prev == Rate
            # assert Conc_citrate - Conc_citrate_Prev == Rate

# Potentially split point for another function

def Plot_Dynamics(Title, Time, Data, Legend):
    X = Time
    Y = Data

    assert len(X) == Y.shape[-1]

    ax = plt.axes(ylim=(0, Y.max() * 1.2))

    ax.set_ylabel('Amount (a.u.)')
    ax.set_xlabel('Time (s)')
    ax.set_title(Title + " over Time")

    for i in range(len(Legend)):
        ax.plot(X, Y[i], label=Legend[i])

    ax.legend(loc='upper left')

    if SaveFilename:
        plt.savefig(SaveFilename)
    else:
        plt.show()

def Indexing(Legend, Query=None):
    SelectedMolecules = list()

    if Query == "TCA":
        # Pick two indices for data
        # Legend = ['oxaloacetate', 'acetyl-CoA', 'H2O', 'citrate', 'CoA', 'AcnA', 'isocitrate', 'Icd', 'NAD+', 'keto-glutarate', 'NADH', 'CO2', 'SucA', 'succinyl-CoA', 'SucD', 'ADP', 'Pi', 'succinate', 'ATP', 'Sdh', 'FAD', 'fumarate', 'FADH2', 'FumA', 'malate', 'Mdh', 'H+']

        SelectedMolecules = ['oxaloacetate', 'acetyl-CoA', 'citrate', 'isocitrate', 'keto-glutarate', 'succinyl-CoA', 'succinate', 'fumarate', 'malate']
        # Subset for convenient visualization, assuming acetyl-CoA level is constant
        SelectedMolecules = ['oxaloacetate', 'citrate', 'isocitrate', 'keto-glutarate', 'succinyl-CoA', 'succinate', 'fumarate', 'malate']
        SelectedMolecules = ['oxaloacetate', 'citrate']

    else:
        if ("S" in Legend) & ("P" in Legend):
            SelectedMolecules = ["S", "P"]
        elif ("S1" in Legend) & ("S2" in Legend):
            SelectedMolecules = ["S1", "S2"]
        elif ("A" in Legend) & ("B" in Legend):
            SelectedMolecules = ["A", "B"]
        elif ("L" in Legend) & ("RL" in Legend) & ("pP" in Legend):
            SelectedMolecules = ["L", "RL", "pP"]
        else:
            SelectedMolecules = Legend

    Idx_SelectedMolecules = list()
    Dict_SelectedMoleculesIdx = dict()
    for SelectedMolecule in SelectedMolecules:
        i = 0
        for Molecule in Legend:
            if Molecule == SelectedMolecule:
                Idx_SelectedMolecules.append(i)
                Dict_SelectedMoleculesIdx[SelectedMolecule] = i
                break
            i += 1

    # PossiblePairs_Legend = [(a, b) for idx_a, a in zip(Idx_SelectedMolecules, SelectedMolecules) for idx_b, b in zip(Idx_SelectedMolecules, SelectedMolecules) if idx_a < idx_b]
    PossiblePairs_Idx = [(idx_a, idx_b) for i, [name_a, idx_a] in enumerate(Dict_SelectedMoleculesIdx.items()) for j, [name_b, idx_b] in enumerate(Dict_SelectedMoleculesIdx.items()) if (i + 1 == j)]
    PossiblePairs_Names = [(name_a, name_b) for i, [name_a, idx_a] in enumerate(Dict_SelectedMoleculesIdx.items()) for j, [name_b, idx_b] in enumerate(Dict_SelectedMoleculesIdx.items()) if (i + 1 == j)]

    return PossiblePairs_Idx, PossiblePairs_Names

def Plot_DynamicsAndPhasePlaneAndRate(Title, Time, Data, Legend, Idx_Pairs, Names_Pairs):
    fig = plt.figure()
    fig.subplots_adjust(wspace=0.2, hspace=0.3)

    N_PossiblePairs = len(Idx_Pairs)
    n = 1
    for i, [a, b] in enumerate(Idx_Pairs):
        X = Data[a]
        Y = Data[b]

        assert Names_Pairs[i][0] == Legend[a]
        assert Names_Pairs[i][1] == Legend[b]

        ax1 = fig.add_subplot(int(N_PossiblePairs), 2, 2 * n - 1)
        ax2 = fig.add_subplot(int(N_PossiblePairs), 2, 2 * n)

        ax1.plot(Time, X, 'r-', label=Legend[a])
        ax1.plot(Time, Y, 'b-', label=Legend[b])
        ax1.set_title(Title + ' Dynamics')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Concentration (a.u.)')
        ax1.legend(loc='upper left')
        ax1.grid()

        # ax = plt.axes(xlim=(0, X.max()), ylim=(0, Y.max()))

        ax2.plot(X, Y, color="blue")
        ax2.set_title(Title + ' Phase plane')
        ax2.set_xlabel(Legend[a])
        ax2.set_ylabel(Legend[b])
        ax2.grid()

        n += 1

    if SaveFilename:
        plt.savefig(SaveFilename)
    else:
        plt.show()


def Plot_Dynamics_Animation(Title, Time, Data, Legend):
    X = Time
    Y = Data

    assert len(X) == Y.shape[-1]

    color = ['red', 'green', 'blue', 'orange']

    fig = plt.figure()
    ax = plt.axes(xlim=(0, 600), ylim=(0, Y.max() * 1.2))

    ax.set_ylabel('Amount (a.u.)')
    ax.set_xlabel('Time (s)')
    ax.set_title(Title + " over Time")

    for i in range(len(Legend)):
        ax.plot(X, Y[i], label=Legend[i])

    ax.legend(loc='upper left')

    lines = plt.plot([], 'b-', markersize=2)
    line = lines[0]

    xanim = []
    yanim = []

    # set data and just scroll through xlims..?
    # line.set_data((x, y))
    offset = 850

    def animate(i):
        index = i + offset
        if index < len(X):
            xanim.append(X[index])
            yanim.append(Y[index])
            ax.set_xlim(min(xanim), index + 10)
            ax.set_ylim(min(yanim) - 10, max(yanim) + 10)

            line.set_data((xanim, yanim))

    ani = FuncAnimation(fig, animate, frames=2000, interval=500)
    plt.show()


def MichaelisMentenEqn(Conc_Enzyme, Conc_Substrate, kcat, KM):
    return (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + KM)



def main():
    Data_Dir = '.'
    Datasets = LoadRawData(Data_Dir)
    for FileName, Dataset in Datasets.items():
        PlotData(Dataset)


if __name__ == '__main__':
    parser = ArgumentParser(description = 'LCC plot')
    parser.add_argument('--save-png',
            dest='save_png',
            type=str,
            help='Plot image file')


    args = parser.parse_args()
    if args.save_png:
        SaveFilename = args.save_png

    main()


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
