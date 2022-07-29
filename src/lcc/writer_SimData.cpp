#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SimData()
{
    std::cout << "Generating simulation state..." << std::endl;

    // write SimData.py
    std::ofstream ofs(Option.SimDataFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "from datetime import datetime" << endl;
    ofs << "import os, sys" << endl;
    ofs << endl;

    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "print('\\n%%%%%% Run " << Option.SimExecutorFile << " to execute simulation. %%%%%%')" << endl;
    ofs << in+ "sys.exit()" << endl;
    ofs << endl;

    ofs << "import numpy as np" << endl;
    ofs << "import csv" << endl;
    ofs << endl;

    // class FDataset
    ofs << "class FDataset:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.Data = None" << endl;
    ofs << endl;

    ofs << in+ "def PrintLegend(self):" << endl;
    ofs << in+ in+ "print(self.Legend)" << endl;
    ofs << endl;

    ofs << in+ "def PrintData(self):" << endl;
    ofs << in+ in+ "print(self.Data)" << endl;
    ofs << endl;

    ofs << in+ "def ExportDataAsList(self):" << endl;
    ofs << in+ in+ "return np.array(self.Data).flatten().tolist()" << endl;
    ofs << endl;

    // class FDataManager
    ofs << "class FDataManager:" << endl;
    ofs << in+ "def __init__(self):" << endl;
    ofs << in+ in+ "# SimOut" << endl;
    ofs << in+ in+ "self.Legend = list()" << endl;
    ofs << in+ in+ "self.DataBuffer = list()" << endl;
    ofs << endl;

    ofs << in+ "def SetLegend(self, InLegend):" << endl;
    ofs << in+ in+ "self.Legend = InLegend" << endl;
    ofs << endl;

    ofs << in+ "def Add(self, InData):" << endl;
    ofs << in+ in+ "self.DataBuffer.append(InData[0])" << endl;
    ofs << endl;

    ofs << in+ "def SaveToTsvFile(self, Legend, Data, InFileName):" << endl;
    ofs << in+ in+ "with open(InFileName, 'w', newline='', encoding='utf-8') as OutFile:" << endl;
    ofs << in+ in+ in+ "TsvWriter = csv.writer(OutFile, delimiter='\\t')" << endl;
    ofs << in+ in+ in+ "if Legend:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(Legend)" << endl;
    ofs << in+ in+ in+ "for Row in Data:" << endl;
    ofs << in+ in+ in+ in+ "TsvWriter.writerow(Row.tolist())" << endl;
    ofs << endl;

    ofs << in+ "def SaveToNpyFile(self, Data, InFileName):" << endl;
    ofs << in+ in+ "np.save('%s.npy' % InFileName, Data)" << endl;
    ofs << endl;

    ofs << in+ "def SaveCountToFile(self, InFileName):" << endl;
    ofs << in+ in+ "self.SaveToTsvFile(self.Legend, self.DataBuffer, InFileName)" << endl;
    ofs << endl;

    std::cout << "  Simulation state has been generated: ";
}
