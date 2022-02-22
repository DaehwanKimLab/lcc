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

bool CheckMolarity(std::string InUnit);
float Prefix2Value(std::string InPrefix);
float Suffix2Value(std::string InSuffix);
float Unit2Value(std::string InUnit);
float Count2Mol(float InCount);
float Mol2Count(float InMol);

} // Namespace Numbers
#endif /* LCC_NUMBER_H */
