#ifndef __LCC_DATA_MANAGER_H__
#define __LCC_DATA_MANAGER_H__

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>

#include "util.h"

class FDataManager {

public:
    void SetLegend(const std::vector<std::string>& InLegend)
    {
        Legend = InLegend;
    }


    void Add(const std::vector<float>& InData) 
    {
        DataBuffer.emplace_back(InData);
    }

    void SaveToFile(const char *InFilename)
    {
        std::ofstream os(InFilename);

        if (!Legend.empty()) {
            std::string Header = Utils::str_join(Legend.begin(), Legend.end(), "\t");
            os << Header << std::endl;
        }

        for (const auto& Row: DataBuffer) {

            std::string tmpstr = Utils::str_join(Row.begin(), Row.end(), "\t");
            os << tmpstr << std::endl;

        }

    }

private:
    std::vector<std::string> Legend;
    std::vector<std::vector<float>> DataBuffer;
};

#endif /* __LCC_DATA_MANAGER_H__ */
