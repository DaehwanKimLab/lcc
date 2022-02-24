#include <iostream>
#include <string>
#include "number.h"

using namespace std;

// avogadro's number
float NA = 6.0221409e+23;

std::map<std::string, float> Prefix2ValueMap {
   { "None", 1 },    // Default 1
   { "Y", 1e24 },    // yotta
   { "Z", 1e21 },    // zetta
   { "E", 1e18 },    // exa
   { "P", 1e15 },    // peta
   { "T", 1e12 },    // tera
   { "G", 1e9 },     // giga
   { "M", 1e6 },     // mega
   { "k", 1e3 },     // kilo
//   { "h", 1e2 },     // hecto
//   { "da", 1e1 },     // deka
   { "d", 1e-1 },    // deci
   { "c", 1e-2 },    // centi
   { "m", 1e-3 },    // milli
   { "u", 1e-6 },    // micro
   { "n", 1e-9 },    // nano
   { "p", 1e-12 },   // pico
   { "f", 1e-15 },   // femto
   { "a", 1e-18 },   // atto
   { "z", 1e-21 },   // zepto
   { "y", 1e-24 },   // yocto
};

std::map<std::string, float> Suffix2ValueMap {
   { "None", 1 },   // Default 1
   { "Other", 1 },   // Default 1
   { "mol", NA },   // mol
   { "M", NA },   // Molarity (treat same as mol, but compiler will handle volume division indication in traversal)
};

namespace Numbers {

bool CheckMolarity(std::string InUnit)
{
    if (InUnit.empty()) {
        return false;
    } else {
        return (InUnit.substr(InUnit.size() - 1) == "M");
    }
}

float Prefix2Value(std::string InPrefix)
{
    return Prefix2ValueMap[InPrefix];
}

float Suffix2Value(std::string InSuffix)
{
    return Suffix2ValueMap[InSuffix];
}

float Unit2Value(std::string InUnit)
{
    auto Unit_Parsed = ParseUnit(InUnit);
    auto Prefix = Unit_Parsed.first;
    auto Suffix = Unit_Parsed.second;

    return (Prefix2Value(Prefix) * Suffix2Value(Suffix));
}

float Unit2ValueWithPrefixOnly(std::string InUnit)
{
    auto Unit_Parsed = ParseUnit(InUnit);
    auto Prefix = Unit_Parsed.first;

    return Prefix2Value(Prefix);
}

std::pair<std::string, std::string> ParseUnit(std::string InUnit)
{
    std::string Prefix;
    std::string Suffix;

    if (InUnit.empty()) {
        Prefix = "None";
        Suffix = "None";
    } else if (InUnit.rfind("mol") != std::string::npos) {
        Prefix = InUnit.substr(0, InUnit.find("mol"));
        Suffix = "mol";
    } else if (CheckMolarity(InUnit)) { // Note: case sensitivity: may need a distinction between meters and Molarity in the future. Cover Mega Molar unit
        Suffix = "M";
        if (InUnit.size() == 1) {
            Prefix = "None";
        } else {
            Prefix = InUnit.substr(0, InUnit.find("M"));
        }
    } else if (InUnit.size() == 1) {
        Prefix = "None";
        Suffix = InUnit;
    } else if (InUnit.size() > 1) {
        Prefix = InUnit[0];
        Suffix = "Other";
    } else {
        Prefix = "None";
        Suffix = "Other";
    }

    std::pair<std::string, std::string> Unit_Parsed(Prefix, Suffix);
    return Unit_Parsed;
}

float Count2Mol(float InCount)
{
    return InCount / NA;
}

float Mol2Count(float InMol)
{
    return InMol * NA;
}

} // Namespace Numbers


// enum Prefix {
//    Y  = 1e24,    // yotta
//    Z  = 1e21,    // zetta
//    E  = 1e18,    // exa
//    P  = 1e15,    // peta
//    T  = 1e12,    // tera
//    G  = 1e9,     // giga
//    M  = 1e6,     // mega
//    k  = 1e3,     // kilo
//    h  = 1e2,     // hecto
//    da = 1e1,     // deka
//    d  = 1e-1,    // deci
//    c  = 1e-2,    // centi
//    m  = 1e-3,    // milli
//    u  = 1e-6,    // micro
//    n  = 1e-9,    // nano
//    p  = 1e-12,   // pico
//    f  = 1e-15,   // femto
//    a  = 1e-18,   // atto
//    z  = 1e-21,   // zepto
//    y  = 1e-24,   // yocto
// }

