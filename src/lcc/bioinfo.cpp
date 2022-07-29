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

std::map<std::string, char> BuildingBlock2Abbr {
        // dNT
        {"dATP", 'A'},  {"dCTP", 'C'},  {"dGTP", 'G'},  {"dTTP", 'T'},
        // NT
        {"ATP", 'A'},   {"CTP", 'C'},   {"GTP", 'G'},   {"TTP", 'U'},
        // AA
        {"ALA", 'A'},   {"ARG", 'R'},   {"ASN", 'D'},   {"ASP", 'N'},   {"CYS", 'C'},   {"GLU", 'E'},   {"GLN", 'Q'},
        {"GLY", 'G'},   {"HIS", 'H'},   {"ILE", 'I'},   {"LEU", 'L'},   {"LYS", 'K'},   {"MET", 'M'},   {"PHE", 'F'},
        {"PRO", 'P'},   {"SER", 'S'},   {"THR", 'T'},   {"TRP", 'W'},   {"TYR", 'Y'},   {"SEC", 'U'},   {"VAL", 'V'},

        // terminator
        {"TER", 'X'}
};

std::map<std::string, std::string> Codon2AA { 
    {"TTT", "PHE"}, {"TTC", "PHE"}, {"TTA", "LEU"}, {"TTG", "LEU"}, 
    {"TCT", "SER"}, {"TCC", "SER"}, {"TCA", "SER"}, {"TCG", "SER"}, 
    {"TAT", "TYR"}, {"TAC", "TYR"}, {"TAA", "TER"}, {"TAG", "TER"}, 
    {"TGT", "CYS"}, {"TGC", "CYS"}, {"TGA", "TER"}, {"TGG", "TRP"}, 
    {"CTT", "LEU"}, {"CTC", "LEU"}, {"CTA", "LEU"}, {"CTG", "LEU"}, 
    {"CCT", "PRO"}, {"CCC", "PRO"}, {"CCA", "PRO"}, {"CCG", "PRO"}, 
    {"CAT", "HIS"}, {"CAC", "HIS"}, {"CAA", "GLN"}, {"CAG", "GLN"}, 
    {"CGT", "ARG"}, {"CGC", "ARG"}, {"CGA", "ARG"}, {"CGG", "ARG"}, 
    {"ATT", "ILE"}, {"ATC", "ILE"}, {"ATA", "ILE"}, {"ATG", "MET"}, 
    {"ACT", "THR"}, {"ACC", "THR"}, {"ACA", "THR"}, {"ACG", "THR"}, 
    {"AAT", "ASN"}, {"AAC", "ASN"}, {"AAA", "LYS"}, {"AAG", "LYS"}, 
    {"AGT", "SER"}, {"AGC", "SER"}, {"AGA", "ARG"}, {"AGG", "ARG"}, 
    {"GTT", "VAL"}, {"GTC", "VAL"}, {"GTA", "VAL"}, {"GTG", "VAL"}, 
    {"GCT", "ALA"}, {"GCC", "ALA"}, {"GCA", "ALA"}, {"GCG", "ALA"}, 
    {"GAT", "ASP"}, {"GAC", "ASP"}, {"GAA", "GLU"}, {"GAG", "GLU"}, 
    {"GGT", "GLY"}, {"GGC", "GLY"}, {"GGA", "GLY"}, {"GGG", "GLY"},
};



//
//std::map<std::string, std::string> BuildingBlock2Abbr {
//        // dNT
//        {"dATP", "A"},  {"dCTP", "C"},  {"dGTP", "G"},  {"dTTP", "T"},
//        // NT
//        {"ATP", "A"},   {"CTP", "C"},   {"GTP", "G"},   {"TTP", "U"},
//        // AA
//        {"ALA", "A"},   {"ARG", "R"},   {"ASN", "D"},   {"ASP", "N"},   {"CYS", "C"},   {"GLU", "E"},   {"GLN", "Q"},
//        {"GLY", "G"},   {"HIS", "H"},   {"ILE", "I"},   {"LEU", "L"},   {"LYS", "K"},   {"MET", "M"},   {"PHE", "F"},
//        {"PRO", "P"},   {"SER", "S"},   {"THR", "T"},   {"TRP", "W"},   {"TYR", "Y"},   {"SEC", "U"},   {"VAL", "V"},
//};
//

//  https://www.bioinformatics.org/sms/iupac.html
//  IUPAC nucleotide code	Base
//  A	Adenine
//  C	Cytosine
//  G	Guanine
//  T (or U)	Thymine (or Uracil)
//  R	A or G
//  Y	C or T
//  S	G or C
//  W	A or T
//  K	G or T
//  M	A or C
//  B	C or G or T
//  D	A or G or T
//  H	A or C or T
//  V	A or C or G
//  N	any base
//  . or -	gap
//

namespace BioInfo {

std::vector<std::string> GetBuildingBlocks(std::string Type)
{
    std::vector<std::string> BuildingBlocks;
    if      ((Type == "dNT") || (Type == "dnt")) { BuildingBlocks = dNT; }
    else if ((Type == "NT") || (Type == "nt"))  { BuildingBlocks = NT; }
    else if ((Type == "AA") || (Type == "aa"))  { BuildingBlocks = AA; }
    else {
        Utils::Assertion(false, "Inappropriate buildingblock type: " + Type);
    }

    return BuildingBlocks;
}

char GetBuildingBlockAbbr(std::string BuildingBlock)
{
    return BuildingBlock2Abbr[BuildingBlock];
}

std::string GetAAFromCodon(std::string Codon)
{
    return Codon2AA[Codon];
}

char GetAAAbbrFromCodon(std::string Codon)
{
    return BuildingBlock2Abbr[Codon2AA[Codon]];
}


} // Namespace BioInfo
