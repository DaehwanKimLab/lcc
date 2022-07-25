import os, sys
from argparse import ArgumentParser
import csv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import numpy as np
import SimFunctions as SimF

SaveFilename = None
NA = 6.0221409e+23
ExclusionString = 'NoPlot'
MolToCount = 'Discrete'
MolToConcentration = 'Conc'

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

    DataCap = 40
    if len(Legend) > DataCap:
        Legend = Legend[:DataCap]
    else:
        DataCap = len(Legend)

    for i, Row in enumerate(Data):
        Data[i] = [float(NumStr) for NumStr in Row[:DataCap]]

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
    # plot_PolymeraseReactions = True
    # plot_Ecoli_NoTCA = True
    # plot_Ecoli_TCA = True
    plot_ShowAll_DynamicsOnly = True
    # plot_ShowAll_DynamicsAndPhasePlain = True

    ### Debugging option ###
    CheckRateExpected = True


    # Detect TCA legend for appropriate display
    if "oxaloacetate" in Legend:
        plot_ShowAll_DynamicsOnly = False
        plot_Ecoli_TCA = True


    Legend_Processed, Data_Processed = ProcessDataToDisplay(Legend, Data_Transposed)

    # Show all
    if plot_ShowAll_DynamicsOnly:
        Cap = 50
        Title = ''
        Data = Data_Processed[1:Cap] / Data_Vol
        Legend = Legend_Processed[1:Cap]

        Plot_Dynamics(Title, Data_Time[0], Data, Legend)
        return 0

    if plot_ShowAll_DynamicsAndPhasePlain:
        Cap = 50
        Title = ''
        Data = Data_Processed[1:Cap] / Data_Vol
        Legend = Legend_Processed[1:Cap]

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

        Data_SimStep = Data_Processed[0:Idx_SimStep]
        Data_Vol = Data_Processed[Idx_SimStep:Idx_Vol]
        Data_PolymeraseCounts = Data_Processed[Idx_Vol:Idx_Polymerase]
        Data_ReactionSubstrateCounts = Data_Processed[Idx_Polymerase:Idx_ReactionSubstrate]
        Data_BuildingBlockCounts = Data_Processed[Idx_ReactionSubstrate:]
        Data_AllSmallMoleculeCounts = Data_Processed[Idx_Polymerase:]

        Legend_SimStep = Legend_Processed[0:Idx_SimStep]
        Legend_Vol = Legend_Processed[Idx_SimStep:Idx_Vol]
        Legend_PolymeraseCounts = Legend_Processed[Idx_Vol:Idx_Polymerase]
        Legend_ReactionSubstrateCounts = Legend_Processed[Idx_Polymerase:Idx_ReactionSubstrate]
        Legend_BuildingBlockCounts = Legend_Processed[Idx_ReactionSubstrate:]
        Legend_AllSmallMoleculeCounts = Legend_Processed[Idx_Polymerase:]

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

        Data_SimStep = Data_Processed[0:Idx_SimStep]
        Data_Vol = Data_Processed[Idx_SimStep:Idx_Vol]
        Data_EnzCounts = Data_Processed[Idx_Vol:Idx_Enzyme]
        Data_SubCounts = Data_Processed[Idx_Enzyme:]

        Legend_SimStep = Legend_Processed[0:Idx_SimStep]
        Legend_Vol = Legend_Processed[Idx_SimStep:Idx_Vol]
        Legend_EnzCounts = Legend_Processed[Idx_Vol:Idx_Enzyme]
        Legend_SubCounts = Legend_Processed[Idx_Enzyme:]

        Title = 'TCA cycle'
        Data = Data_SubCounts
        Legend = Legend_SubCounts

    elif plot_TCACycle_Reduced:
        Molecules = ['SimStep', 'Vol', 'GltA', 'oxaloacetate', 'acetyl-CoA', 'H2O', 'citrate', 'CoA']

        Idx_SimStep = 1
        Idx_Vol = 1 + Idx_SimStep
        Idx_Enzyme = 1 + Idx_Vol

        Data_SimStep = Data_Processed[0:Idx_SimStep]
        Data_Vol = Data_Processed[Idx_SimStep:Idx_Vol]
        Data_EnzCounts = Data_Processed[Idx_Vol:Idx_Enzyme]
        Data_SubCounts = Data_Processed[Idx_Enzyme:]

        Legend_SimStep = Legend_Processed[0:Idx_SimStep]
        Legend_Vol = Legend_Processed[Idx_SimStep:Idx_Vol]
        Legend_EnzCounts = Legend_Processed[Idx_Vol:Idx_Enzyme]
        Legend_SubCounts = Legend_Processed[Idx_Enzyme:]

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

        Data = Data_Processed[Target:Target+N_Species]
        Legend = Legend_Processed[Target:Target+N_Species]

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

        Data = Data_Processed[Target:Target+N_Species]
        Legend = Legend_Processed[Target:Target+N_Species]

    #
    # # All TCA Cases
    # Idx_Pairs, Names_Pairs = Indexing(Legend, Query="TCA")
    #
    # if CheckRateExpected:
    #     RateCheck(Data, Legend, Idx_Pairs, Names_Pairs)
    #
    # Plot_DynamicsAndPhasePlaneAndRate(Title, Data_Time[0], Data, Legend, Idx_Pairs, Names_Pairs)


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

    # Data display filtering
    IdxToDelete = list()
    for i in range(len(Legend)):
        if Legend[i] == "Pseudo":
            IdxToDelete.append(i)
        if Legend[i] == "H2O" or Legend[i] == "CO2":
            IdxToDelete.append(i)

    if IdxToDelete:
        Legend = [j for i, j in enumerate(Legend) if i not in IdxToDelete]
        Y = np.delete(Y, np.array(IdxToDelete), axis=0)

    unit_txt = 'a.u.'

    array_unit = {
        'nM' : 1e-9,
        'uM' : 1e-6,
        'mM' : 1e-3
    }

    for utxt, uval in array_unit.items():
        if np.any(Y > uval):
            print('Unit has been set to', utxt)
            unit_txt = utxt

    if unit_txt in array_unit:
        Y = Y / NA / array_unit[unit_txt]
        print('Final unit:', unit_txt)

    ax = plt.axes(xlim=(0, X.max()), ylim=(0, Y.max() * 1.2))

    ax.set_ylabel('Amount (' + unit_txt + ')')
    ax.set_xlabel('Time (s)')
    ax.set_title(Title + " over Time")

    for i in range(len(Legend)):
        line, = ax.plot(X, Y[i], label=Legend[i])
        ax.legend(loc='upper left')
        if len(Legend) > 4:
            ax.text(X[-1] * 1.01, Y[i][-1], Legend[i] + ": {}".format(Y[i][-1]), va="center", color=line.get_color())

    if SaveFilename:
        plt.savefig(SaveFilename)
    else:
        plt.show()

def ProcessDataToDisplay(Legend, Data):
    Data_Processed = list()
    Legend_Processed = list()
    Size = len(Data[0])
    for i, (Data_Row, Legend_Element) in enumerate(zip(Data, Legend)):
        if Legend_Element == 'Vol':
            continue
        if np.array_equal(Data_Row, np.full_like(Data_Row, Data_Row[0])):
            continue
        if ExclusionString in Legend_Element:
            continue
        if MolToCount in Legend_Element:
            Data_Row *= NA
        if MolToConcentration in Legend_Element:
            Data_Row /= NA
        Data_Processed.append(Data_Row)
        Legend_Processed.append(Legend_Element)

    return Legend_Processed, np.array(Data_Processed)

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


# def Plot_Dynamics_Animation(Title, Time, Data, Legend):
#     X = Time
#     Y = Data
#
#     assert len(X) == Y.shape[-1]
#
#     color = ['red', 'green', 'blue', 'orange']
#
#     fig = plt.figure()
#     ax = plt.axes(xlim=(0, 600), ylim=(0, Y.max() * 1.2))
#
#     ax.set_ylabel('Amount (a.u.)')
#     ax.set_xlabel('Time (s)')
#     ax.set_title(Title + " over Time")
#
#     for i in range(len(Legend)):
#         ax.plot(X, Y[i], label=Legend[i])
#
#     ax.legend(loc='upper left')
#
#     ani = FuncAnimation(fig, animate, frames=2000, interval=500)
#     plt.show()
#
#
#
# lines = plt.plot([], 'b-', markersize=2)
# line = lines[0]
#
# xanim = []
# yanim = []
#
# # set data and just scroll through xlims..?
# # line.set_data((x, y))
# offset = 850
#
#
# def animate(i):
#     index = i + offset
#     if index < len(x):
#         xanim.append(x[index])
#         yanim.append(y[index])
#         ax.set_xlim(min(xanim), index + 10)
#         ax.set_ylim(min(yanim) - 10, max(yanim) + 10)
#
#         line.set_data((xanim, yanim))
#
#
# ani = FuncAnimation(fig, animate, frames=2000, interval=500)
# plt.show()
#

def MichaelisMentenEqn(Conc_Enzyme, Conc_Substrate, kcat, KM):
    return (kcat * Conc_Enzyme * Conc_Substrate) / (Conc_Substrate + KM)



def main(SaveToFile = None):
    global SaveFilename
    if SaveToFile is not None:
        SaveFilename = SaveToFile

    Data_Dir = '.'
    Datasets = LoadRawData(Data_Dir)
    for FileName, Dataset in Datasets.items():
        PlotData(Dataset)




def Plot3D(Nodes, dim=None, distance=None, volume=None, shape='cuboid'):
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    X = Nodes[:, 0]
    Y = Nodes[:, 1]
    Z = Nodes[:, 2]
    ax.plot3D(X, Y, Z, 'red')

    if np.any(dim):
        ax.set_box_aspect(dim)
    #
    # ax.set_xlim3d(0, X.shape[0])
    # ax.set_ylim3d(0, Y.shape[0])
    # ax.set_zlim3d(0, Z.shape[0])

    # Data for three-dimensional scattered points
    # ax.scatter3D(X, Y, Z, c='blue')
    # ax.scatter3D(X, Y, Z, c=Z, cmap='Greens')

    VolTxt = ''

    if np.any(dim):
        VolTxt = 'Volume: '
        Vol = 0
        if shape == 'cuboid':
            Vol = np.prod(np.array(dim)) / 1e9
        elif shape == 'ellipsoid':
            Vol = (4 / 3) * np.pi * np.prod(np.array(dim) / 2) / 1e9
        elif shape == 'cylinder':
            Vol = np.pi * np.prod(np.array([dim[0] / 2, dim[1], dim[2] / 2])) / 1e9
        VolTxt += str(Vol) + 'um^3\n '

    ax.set_title(VolTxt + 'Distance Covered: {:.3f} nm\n # of Nodes (retained): {}'.format(distance, Nodes.shape[0]))

    plt.show()


if __name__ == '__main__':
    parser = ArgumentParser(description = 'LCC plot')
    parser.add_argument('--save-fig',
            dest='save_fname',
            type=str,
            help='Save figure to file')

    args = parser.parse_args()

    main(args.save_fname)


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
