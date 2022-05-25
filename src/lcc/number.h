#ifndef LCC_NUMBER_H
#define LCC_NUMBER_H

#include <string>
#include <vector>
#include <map>

namespace Numbers {

// SI base units
// length	meter	m
// mass	kilogram	kg
// time	second	s
// electric current	ampere	A
// amount of substance	mole	mol
// luminous intensity candela	cd

float GetAvogadro();
float MultiplyByAvogadro(float InFloat);
std::string GetAvogadroStr();
float RandomNumber(float Min, float Max);
float GetFloatDefault();
float GetIntDefault();
bool CheckMolarity(std::string Unit);
float Prefix2Value(std::string Prefix);
float Suffix2Value(std::string Suffix);
float Unit2Value(std::string Unit);
float Unit2ValueWithPrefixOnly(std::string Unit);
std::pair<std::string, std::string> ParseUnit(std::string Unit);
float Conversion_bp2nm(int BP);
float Count2Mol(float Count);
float Mol2Count(float Mol);

} // Namespace Numbers
#endif /* LCC_NUMBER_H */
