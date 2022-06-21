import csv
import numpy as np

from pypdb import *

# for advanced search
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.search_client import QueryGroup, LogicalOperator
from pypdb.clients.search.operators import text_operators


def LoadTSVDatabase(db_fname):
    db = None
    with open(db_fname) as fp:
        csv_reader = csv.reader(fp, delimiter='\t')
        list_of_rows = list(csv_reader)
        db = list_of_rows[1:]
    return db


def OpenTSVDatabase(db_fname):
    db = LoadTSVDatabase(db_fname)
    Database_Gene = ParseGenes(db)

    print('Database constructed')

    return Database_Gene


def SaveDictToTSV(Database, KeysToSave, InFileName):
    Legend = KeysToSave

    with open(InFileName + '.tsv', 'w', newline='', encoding='utf-8') as OutFile:
        TsvWriter = csv.writer(OutFile, delimiter='\t')
        if Legend:
            TsvWriter.writerow(Legend)

        for i in range(len(Database[KeysToSave[0]])):
            Row = list()
            for Key in KeysToSave:
                Row.append(Database[Key][i])
            TsvWriter.writerow(Row)





def ParseGenes(db_genes):
    db = dict()
    NUniq_Genes = len(db_genes)
    db['Symbol'] = list()
    db['Length'] = np.zeros(NUniq_Genes)
    db['Coord'] = np.zeros(NUniq_Genes)
    db['Dir'] = np.zeros(NUniq_Genes)
    db['Seq'] = list()
    db['PDBID'] = list()
    db['GeneToPDBID'] = list()

    Dir = dict()
    Dir['+'] = 1
    Dir['-'] = -1

    for i, Value in enumerate(db_genes):
        Length, Name, Seq, RNAID, Coordinate, Direction, Symbol, Type, GeneID, MonomerID = Value
        db['Symbol'].append(Symbol)
        # db['Length'][i] = (len(Seq))
        # db['Coord'][i] = int(Coordinate)
        # db['Dir'][i] = Dir[Direction]
        # db['Seq'].append(Seq)

        PDBID = SearchPDB_OldAPI(Symbol)
        db['PDBID'].append(PDBID)
        # db['GeneToPDBID'].append((Symbol, PDBID))

        if i % 100 == 0:
            print('# of genes scanned: ', i)

        # DL: Debugging Purposes
        # if i == 3:
        #     break

    print('Total # of genes parsed: ', len(db['Symbol']))

    return db


def SearchPDB_OldAPI(Symbol):
    found_pdbs = Query(Symbol).search()

    if found_pdbs:
        PDBID = found_pdbs[0]
        # print(Symbol, found_pdbs)
    else:
        PDBID = ''
        # print(Symbol, [PDBID])

    # print(PDBID)

    return PDBID


def SimpleSearchPDB(Symbol):
    search_operator = text_operators.DefaultOperator(value=Symbol)
    return_type = ReturnType.ENTRY

    results = perform_search(search_operator, return_type)

    print(Symbol, results[:10])

    PDBID = results[0]

    return PDBID


def AdvSearchPDB(Symbol, Organism):
    search_operator = text_operators.DefaultOperator(value=Symbol)

    is_organism_operator = text_operators.ExactMatchOperator(
        value=Organism,
        attribute="rcsb_entity_source_organism.taxonomy_lineage.name")

    search_and_is_organism_operator = QueryGroup(
        queries=[search_operator, is_organism_operator],
        logical_operator=LogicalOperator.AND
    )

    return_type = ReturnType.ENTRY

    results = perform_search_with_graph(
        query_object=search_and_is_organism_operator,
        return_type=return_type)

    print(Symbol, results[:10])  # Huzzah

    PDBID = results[0]

    return PDBID


def main():
    DatabaseFileName = r'./Database/genes.tsv'
    Database = OpenTSVDatabase(DatabaseFileName)

    Genes = Database['Symbol']
    PDBIDs = Database['PDBID']

    KeysToSave = ['Symbol', 'PDBID']

    KeyToSave = ['GeneToPDBID']
    SaveFileName = 'GenePDBID'

    SaveDictToTSV(Database, KeysToSave, SaveFileName)


if __name__ == '__main__':
    main()

