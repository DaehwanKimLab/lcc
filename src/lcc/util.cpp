#ifdef _MSC_VER
#include <Windows.h>
#include <fileapi.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>


#include "util.h"

using namespace std;

namespace Utils {

bool CreatePath(const char *Path)
{
#ifdef _MSC_VER
	return CreateDirectoryA(Path, NULL);
#else
    return mkdir(Path, 0775) == 0 || (errno == EEXIST);
#endif
}

bool CreatePaths(const char *Path)
{
    vector<string> Fields;

#ifdef _MSC_VER
	const char* DirectorySep = "\\";
#else
	const char* DirectorySep = "/";
#endif

    tokenize(Path, DirectorySep, Fields);


    if (Fields.size() == 0) {
        return true;
    }

    string NewPath;

    for(const auto& p : Fields) {
        NewPath += p + DirectorySep;

        if (!CreatePath(NewPath.c_str())) {
            return false;
        }
    }

    return true;
}

std::string SciFloat2Str (float Float)
{
   std::stringstream ss;
   ss << std::scientific << Float;
   return ss.str();
}

}
