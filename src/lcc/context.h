#ifndef LCC_CONTEXT_H
#define LCC_CONTEXT_H

#include <string>
#include <vector>
#include <map>

#include "option.h"

typedef std::map<std::string, std::string> FTableRecord;

class FTable {
public:
	std::vector<FTableRecord> Records;

	void LoadFromTSV(const char *Filename);
	void Dump();

};


class FCompilerContext {
public:

    void Init(const FOption& InOption);

    FTable GeneTable;
};

#endif /* LCC_CONTEXT_H */
