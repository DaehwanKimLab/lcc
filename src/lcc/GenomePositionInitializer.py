# from mpl_toolkits import mplot3d

# %matplotlib inline
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
from sqlite3 import Error
import SimFunctions as SimF

# Temporary
import csv

'''
For initializing Genome location in a 3D cell container
'''


def DetermineDistance_3D(P1, P2):
    return np.sqrt(np.sum((P1 - P2) ** 2, axis=1))


def Convert_bp2nm(BP):
    return (BP / 10) * 3.4


def Convert_nm2bp(NM):
    return (NM / 3.4) * 10


def TrimToShape(Nodes, Dim, Shape):
    if Shape == 'cuboid':
        return Nodes
    elif Shape == 'ellipsoid':
        # ellipsoid equation: x**2/a**2 + y**2/b**2 + z**2/c**2 = 1, where a, b, c are the principal semiaxes.
        Semiaxes = np.array([Dim[0], Dim[1], Dim[2]]) / 2
        IdxToRetain = np.where(np.sum(((Nodes - Semiaxes) ** 2) / (Semiaxes ** 2), axis=1) < 1)
        return Nodes[IdxToRetain]
    elif Shape == 'cylinder':
        # 2 π r^2 are the principal semiaxes.
        Radius = Dim[0] / 2
        Ref_XZ = np.array([Radius, Radius])   # may just use broadcasting
        DistancesFromRef = DetermineDistance_3D(Ref_XZ, np.stack([Nodes[:, 0], Nodes[:, 2]], axis=1))
        IdxToRetain = np.where(DistancesFromRef < Radius)
        return Nodes[IdxToRetain]
    else:
        print('Unrecognizable Shape for Trimming: ', Shape)
        print('Available Trimming Shapes: cuboid, ellipsoid, cylinder')
        sys.exit(0)


def GetNodesAndDistances(Dim, N_Nodes, shape='cuboid'):
    assert (isinstance(N_Nodes, int)) & (
                N_Nodes > 1), 'ERROR: N_Nodes (Input: %s) must be an integer greater than 1.' % N_Nodes

    XYZ = np.array([Dim])

    # Generate sets of three random numbers
    Nodes_Original = np.random.rand(N_Nodes, 3) * XYZ

    # TODO:Potential node trimming step for shape control
    Nodes_Trimmed = TrimToShape(Nodes_Original, Dim, shape)

    # Arrange Nodes by finding the nearest neighbor
    Node_Reference = Nodes_Trimmed[0]
    Nodes_Modified = Nodes_Trimmed[1:]  # Nodes excluding the reference node

    Nodes_Arranged = np.zeros_like(Nodes_Trimmed)
    Nodes_Arranged[0] = Node_Reference

    Distances = np.zeros(Nodes_Trimmed.shape[0])

    # print('Prev', Nodes_Modified)
    # print('Node:', Node_Reference)

    for i in range(Nodes_Modified.shape[0]):
        Distance = DetermineDistance_3D(Node_Reference, Nodes_Modified)
        Idx_NearestNeighbor = Distance.argmin()
        Distances[i] = Distance[Idx_NearestNeighbor]
        Node_Reference = Nodes_Modified[Idx_NearestNeighbor]
        Nodes_Arranged[i + 1] = Node_Reference
        Nodes_Modified = np.delete(Nodes_Modified, Idx_NearestNeighbor, axis=0)
        # print(i)
        # print('Node    :', Node_Reference)
        # print('Rest    :', Nodes_Modified)
        # print('Arranged:', Nodes_Arranged)

    return Nodes_Arranged, Distances


def Plot3D(Nodes, Distance, dim=None, shape='cuboid'):
    ax = plt.axes(projection='3d')

    # Data for a three-dimensional line
    X = Nodes[:, 0]
    Y = Nodes[:, 1]
    Z = Nodes[:, 2]
    ax.plot3D(X, Y, Z, 'red')

    if dim:
        ax.set_box_aspect(dim)
    #
    # ax.set_xlim3d(0, X.shape[0])
    # ax.set_ylim3d(0, Y.shape[0])
    # ax.set_zlim3d(0, Z.shape[0])

    # Data for three-dimensional scattered points
    # ax.scatter3D(X, Y, Z, c='blue')
    # ax.scatter3D(X, Y, Z, c=Z, cmap='Greens')

    VolTxt = ''

    if dim:
        VolTxt = 'Volume: '
        Vol = 0
        if shape == 'cuboid':
            Vol = np.prod(np.array(dim)) / 1e9
        elif shape == 'ellipsoid':
            Vol = (4 / 3) * np.pi * np.prod(np.array(dim) / 2) / 1e9
        elif shape == 'cylinder':
            Vol = np.pi * np.prod(np.array([dim[0] / 2, dim[1], dim[2] / 2])) / 1e9
        VolTxt += str(Vol) + 'um^3\n '

    ax.set_title(VolTxt + 'Distance Covered: {:.3f} nm\n # of Nodes (retained): {}'.format(Distance, Nodes.shape[0]))

    plt.show()


# def SetUpDatabase(Database_fname, type='tsv'):
#     Database = None
#
#     if type == 'tsv':
#         Database = OpenTSVDatabase(Database_fname)
#         return Database
#
#     elif type == 'sql':
#         Database = OpenSQLDatabase(Database_fname)
#         return Database
#
#     else:
#         print("Error: Database type does not exist: %s" % type, file=sys.stderr)
#         sys.exit(0)


# Temporary TSV Database


def LoadTSVDatabase(db_fname):
    db = None

    with open(db_fname) as fp:
        csv_reader = csv.reader(fp, delimiter='\t')
        list_of_rows = list(csv_reader)
        db = list_of_rows[1:]

    return db


def ParseGenes(db_genes):
    db = dict()
    NUniq_Genes = len(db_genes)
    db['Symbol'] = list()
    db['Length'] = np.zeros(NUniq_Genes)
    db['Coord'] = np.zeros(NUniq_Genes)
    db['Dir'] = np.zeros(NUniq_Genes)
    db['Seq'] = list()

    Dir = dict()
    Dir['+'] = 1
    Dir['-'] = -1

    for i, Value in enumerate(db_genes):
        Length, Name, Seq, RNAID, Coordinate, Direction, Symbol, Type, GeneID, MonomerID = Value
        db['Symbol'].append(Symbol)
        db['Length'][i] = (len(Seq))
        db['Coord'][i] = int(Coordinate)
        db['Dir'][i] = Dir[Direction]
        db['Seq'].append(Seq)

    return db


def OpenTSVDatabase(db_fname):
    db = None

    # Load
    db = LoadTSVDatabase(db_fname)

    # Parse
    Database_Gene = ParseGenes(db)

    return Database_Gene


# From EcoliLanguage/Scripts/CreateSQLiteDatabase.py
def OpenSQLDatabase(db_fname):
    db = None

    try:
        db = sqlite3.connect(db_fname)
        return db
    except sqlite3.Error as e:
        print("Error: {}".format(e), file=sys.stderr)
        sys.exit(0)


def OpenDatabase(db_fname):
    if '.tsv' in db_fname:
        return OpenTSVDatabase(db_fname)
    else:
        return OpenSQLDatabase(db_fname)


def ExecuteSQLQuery(db, sql_query):
    c = db.cursor()
    c.execute(sql_query)
    return c.fetchall()


def GetXYZForGenomePositionsInBP(Positions_bp, Nodes, Distances):
    Positions_nm = (Positions_bp / 10) * 3.4   # unit conversion from bp to nm
    NodePositionOnGenome = np.reshape(np.cumsum(np.roll(Distances, shift=1)), [1, -1])

    # # For Debugging
    # Positions_nm = Positions_nm[0:4]
    # NodePositionOnGenome *= 2000   # for debugging

    Positions_Bin = Positions_nm < NodePositionOnGenome
    Positions_Bin_Cumsum = np.cumsum(Positions_Bin, axis=1)
    Positions_Idx = np.where(Positions_Bin_Cumsum == 1)[1]

    Point_1 = NodePositionOnGenome[:, Positions_Idx - 1]
    Point_2 = NodePositionOnGenome[:, Positions_Idx ]

    PercentLengthOfGeneStartBetweenPoints = ((Positions_nm.transpose() - Point_1) / (Point_2 - Point_1)).transpose()

    XYZ_1 = Nodes[Positions_Idx - 1]
    XYZ_2 = Nodes[Positions_Idx]

    Positions_XYZ = XYZ_1 + ((XYZ_2 - XYZ_1) * PercentLengthOfGeneStartBetweenPoints)
    assert np.count_nonzero(Positions_XYZ < 0) == 0, 'ERROR: XYZ coordinates cannot be less than zero'

    return Positions_XYZ


def main():   # add verbose

    # In nano meter scale assumption

    # Ecoli info

    # 0.6–0.7 μm3 in volume, cylinder 1.0-2.0 um long, with radius 0.5 um.
    Dim_X = 800
    Dim_Y = 1200    # 1000 ~ 2000 nm
    Dim_Z = 800
    Dim = (Dim_X, Dim_Y, Dim_Z)
    N_Nodes = 50000

    ChrSize_bp = 4641652
    ChrSize_nm = Convert_bp2nm(ChrSize_bp)    # 1,578,161.68 nm

    # Shape = 'ellipsoid'
    Shape = 'cylinder'

    if Shape == 'cylinder':
        assert Dim_X == Dim_Z, 'cylinder shape must have same X and Z dimensions as a radius'

    # Genome Position Location
    Nodes, Distances = GetNodesAndDistances(Dim, N_Nodes, shape=Shape)

    # Visualization of the result
    Plot3D(Nodes, np.sum(Distances), Dim, shape=Shape)


    # Set database filename
    # DatabaseFileName = ''
    DatabaseFileName = r'./Database/genes.tsv'

    # Open Database for annotation (temporary options: 'tsv')
    Database = OpenDatabase(DatabaseFileName)

    Gene_Start_bp = np.reshape(Database['Coord'], [-1, 1])
    Gene_End_bp = np.reshape(Database['Coord'] + Database['Length'] * Database['Dir'], [-1, 1])

    # Output

    # ChrSize = ChrSize_bp
    # ChrSeq
    # ChrNodes = Nodes

    Gene_Names = Database['Symbol']
    Gene_Start_XYZ = GetXYZForGenomePositionsInBP(Gene_Start_bp, Nodes, Distances)
    Gene_End_XYZ = GetXYZForGenomePositionsInBP(Gene_End_bp, Nodes, Distances)


    # # Get Annotation
    # Query = '''
    # SELECT * FROM Gene;
    #
    # '''
    # Result = ExecuteSQLQuery(Database, Query)
    # print(Result)
    #
    # # Align Annotation
    #

if __name__ == '__main__':
    main()



