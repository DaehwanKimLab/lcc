#ifndef LCC_STATE_H
#define LCC_STATE_H

#include <string>
#include <vector>
#include <map>


#include "datamanager.h"
// #include "dataset.h"

extern FDataManager DataManager;

class FState {
public:

// private:
    std::map<std::string, int> Enz2Count; // The count info may be merged but in the future sim.py
    std::map<std::string, int> Sub2Count;
    float Vol;
    float SimStep = 0; // temporary

    void SetVol(const int& InVol);
    void SetEnzCount(const std::string& Name, const int& Count);
    void SetSubCount(const std::string& Name, const int& Count);
    int GetEnzCount(const std::string& Name);
    int GetSubCount(const std::string& Name);
    float GetEnzConc(const std::string& Name);
    float GetSubConc(const std::string& Name);
    void PrintState();
    std::vector<std::string> ExportLegend();
    std::vector<float> ExportState();
};
#endif /* LCC_STATE_H */
