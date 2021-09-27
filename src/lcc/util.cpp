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
    return mkdir(Path, 0775) == 0 || (errno == EEXIST);
}

bool CreatePaths(const char *Path)
{
    vector<string> Fields;

    tokenize(Path, "/", Fields);


    if (Fields.size() == 0) {
        return true;
    }

    string NewPath;

    for(const auto& p : Fields) {
        NewPath += p + "/";

        if (!CreatePath(NewPath.c_str())) {
            return false;
        }
    }

    return true;
}

}
