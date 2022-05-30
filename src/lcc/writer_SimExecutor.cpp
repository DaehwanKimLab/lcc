#include "writer.h"
#include <algorithm>

using namespace std;

void FWriter::SimExecutor()
{
    std::cout << std::endl << "Generating simulation executor..." << std::endl;

    // write SimModule.py
    std::ofstream ofs(Option.SimExecutorFile.c_str());
    std::string endl = "\n";

    // IMPORT
    ofs << "import os, sys" << endl;
    ofs << "from argparse import ArgumentParser" << endl;
    ofs << endl;
    // ofs << "import SimIdx as idx" << endl;
    if (Context.LocationList.empty())   { ofs << "import SimModule" << endl;
                                          ofs << "import plot" << endl;      }

    else                                { ofs << "import SimVis2D" << endl;  }
                                        { ofs << "import SimServer" << endl;  }
    ofs << endl;

    // BODY
    ofs << "def main():   # add verbose" << endl;
    // Execute appropriate simulation
    if (Context.LocationList.empty())   { ofs << in+ "SimModule.main()" << endl;
                                          ofs << in+ "plot.main(args.save_fname)" << endl; }

    else                                { ofs << in+ "SimServer.main()" << endl; }
                                        { ofs << in+ "SimVis2D.main()" << endl; }
    ofs << endl;

    // MAIN
    ofs << "if __name__ == '__main__':" << endl;
    ofs << in+ "parser = ArgumentParser()" << endl;
    ofs << in+ "parser.add_argument('--save-fig', dest='save_fname', type=str, help='Save figure to file')" << endl;
    ofs << in+ "args = parser.parse_args()" << endl;
    ofs << in+ "main()" << endl;
    ofs << endl;

    std::cout << "  Simulation executor has been generated: ";
}
