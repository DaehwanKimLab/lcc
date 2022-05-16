#ifndef LCC_BIOINFO_H
#define LCC_BIOINFO_H

#include <string>
#include <vector>
#include <map>

namespace BioInfo {

std::vector<std::string> GetBuildingBlocks(std::string Type);
char GetBuildingBlockAbbr(std::string BuildingBlock);
std::string GetAAFromCodon(std::string Codon);

} // Namespace BioInfo
#endif /* LCC_BIOINFO_H */
