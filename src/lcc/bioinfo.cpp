#include <iostream>
#include <string>
#include "bioinfo.h"
#include "util.h"
#include "number.h"
#include <random>

using namespace std;

std::vector<std::string> dNT = {"dATP", "dCTP", "dGTP", "dTTP"};
std::vector<std::string> NT = {"ATP", "CTP", "GTP", "UTP"};
std::vector<std::string> AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLT", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "SEL", "VAL"};

std::map<std::string, std::string> BuildingBlock2Abbr {
        // dNT
        {"dATP", "A"},  {"dCTP", "C"},  {"dGTP", "G"},  {"dTTP", "T"},
        // NT
        {"ATP", "A"},   {"CTP", "C"},   {"GTP", "G"},   {"TTP", "U"},
        // AA
        {"ALA", "A"},   {"ARG", "R"},   {"ASN", "D"},   {"ASP", "N"},   {"CYS", "C"},   {"GLU", "E"},   {"GLN", "Q"},
        {"GLY", "G"},   {"HIS", "H"},   {"ILE", "I"},   {"LEU", "L"},   {"LYS", "K"},   {"MET", "M"},   {"PHE", "F"},
        {"PRO", "P"},   {"SER", "S"},   {"THR", "T"},   {"TRP", "W"},   {"TYR", "Y"},   {"SEC", "U"},   {"VAL", "V"},
};


namespace BioInfo {

std::vector<std::string> GetBuildingBlocks(std::string Type)
{
    std::vector<std::string> BuildingBlocks;
    if      ((Type == "dNT") || (Type == "dnt")) { BuildingBlocks = dNT; }
    else if ((Type == "NT") || (Type == "nt"))  { BuildingBlocks = NT; }
    else if ((Type == "AA") || (Type == "aa"))  { BuildingBlocks = AA; }
    else                    { Utils::Assertion(false, "Inappropriate buildingblock type: " + Type); }

    return BuildingBlocks;
}

std::string GetBuildingBlockAbbr(std::string BuildingBlock)
{
    return BuildingBlock2Abbr[BuildingBlock];
}


} // Namespace BioInfo
