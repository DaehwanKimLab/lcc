#ifdef _MSC_VER
#include <Windows.h>
#include <fileapi.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
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

std::string SciFloat2Str(float Float)
{
   std::stringstream ss;
   ss << std::scientific << Float;
   return ss.str();
}

std::string JoinStr2Str(std::vector<std::string> StringList)
{
    std::string JoinedStr;
    for (auto& Str : StringList){
        JoinedStr += "'" + Str + "'" + ", ";
    }
    return JoinedStr;
}

std::string JoinInt2Str(std::vector<int> IntList)
{
    std::string JoinedStr;
    for (auto& Int : IntList){
        JoinedStr += std::to_string(Int) + ", ";
    }
    return JoinedStr;
}

std::string JoinInt2Str_Idx(std::vector<int> IntList)
{
    std::string JoinedStr;
    for (auto& Int : IntList){
        JoinedStr += std::to_string(Int) + ", ";
    }
    return JoinedStr;
}

std::string JoinFloat2Str(std::vector<float> FloatList)
{
    std::string JoinedStr;
    for (auto& Float : FloatList){
        JoinedStr += Utils::SciFloat2Str(Float) + ", ";
    }
    return JoinedStr;
}

std::string Matrix2Str(std::vector<std::vector<int>> Matrix)
{
    std::string MatrixStr;
    int N_Rows;
    int N_Columns;

    N_Rows = Matrix.size();
    for (auto& Row : Matrix) {
        MatrixStr += "[";
        N_Columns = Row.size();
        for (auto& Item : Row) {
            MatrixStr += std::to_string(Item) + ", ";
        }
        MatrixStr += "], ";
    }
    return MatrixStr;
}

void Assertion(bool Bool, std::string ErrorMessage)
{
    if (!Bool) {
        std::cout << "%%%%%%  " << ErrorMessage << "  %%%%%%" << std::endl;
        assert(Bool);
        }
}

}
