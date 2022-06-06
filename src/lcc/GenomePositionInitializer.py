# from mpl_toolkits import mplot3d

# %matplotlib inline
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import sqlite3
from sqlite3 import Error
import SimFunctions as SimF
import plot

# Temporary
import csv

'''
For initializing Genome location in a 3D cell container
'''

def Convert_bp2nm(BP):
    return (BP / 10) * 3.4


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


def main():   # add verbose

    # In nano meter scale assumption

    print('################# TEST GenomePositionInitializer #################')

    # Ecoli info

    # 0.6–0.7 μm3 in volume, cylinder 1.0-2.0 um long, with radius 0.5 um.
    Dim_X = 800
    Dim_Y = 1200    # 1000 ~ 2000 nm
    Dim_Z = 800
    Dim = (Dim_X, Dim_Y, Dim_Z)

    ChrSize_bp = 4641652
    ChrLen_nm = Convert_bp2nm(ChrSize_bp)    # 1,578,161.68 nm

    # Shape = 'ellipsoid'
    Shape = 'cylinder'

    if Shape == 'cylinder':
        assert Dim_X == Dim_Z, 'cylinder shape must have same X and Z dimensions as a radius'

    print('Dimensions: ', Dim)
    print('ChrSize_bp: ', ChrSize_bp)
    print('ChrLen_nm:  ', ChrLen_nm)
    print('Shape:      ', Shape)

    N_Nodes = 150000
    print('N_Nodes:    ', N_Nodes)

    # Genome Position Location
    Nodes, Distances = SimF.GetNodesAndDistances(Dim, ChrLen_nm, shape=Shape, n_nodes=N_Nodes)

    # Save if sufficient for the genome length
    TotalDistance = np.sum(Distances)
    if TotalDistance > ChrLen_nm:
        np.save('./Database/EcoliGenomeNodes.npy', Nodes)
        np.save('./Database/EcoliGenomeDistances.npy', Distances)

    # Set database filename
    # DatabaseFileName = ''
    DatabaseFileName = r'./Database/genes.tsv'

    # Open Database for annotation (temporary options: 'tsv')
    Database = OpenTSVDatabase(DatabaseFileName)

    Gene_Start_bp = np.reshape(Database['Coord'], [-1, 1])
    Gene_End_bp = np.reshape(Database['Coord'] + Database['Length'] * Database['Dir'], [-1, 1])

    # Output

    # ChrSize = ChrSize_bp
    # ChrSeq
    # ChrNodes = Nodes

    Gene_Names = Database['Symbol']
    Gene_Start_XYZ = SimF.GetXYZForGenomePositionsInBP(Gene_Start_bp, Nodes, Distances)
    Gene_End_XYZ = SimF.GetXYZForGenomePositionsInBP(Gene_End_bp, Nodes, Distances)

    print('Gene_Names', Gene_Names)
    print('Gene_Start_XYZ', Gene_Start_XYZ)
    print('Gene_End_XYZ', Gene_End_XYZ)

    # Visualization of the result
    plot.Plot3D(Nodes, dim=Dim, distance=TotalDistance, shape=Shape)

    print('')


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



