#ifndef LCC_UTIL_H
#define LCC_UTIL_H
#include <iterator>
#include <sstream>
/**
 * Copied from Bowtie2's implementation
 * Split string s according to given delimiters.  Mostly borrowed
 * from C++ Programming HOWTO 7.3.
 */

namespace Utils {
template<typename T>
static inline void tokenize(
        const std::string& s,
        const std::string& delims,
        T& ss,
        size_t max = -1)
{
    //string::size_type lastPos = s.find_first_not_of(delims, 0);
    std::string::size_type lastPos = 0;
    std::string::size_type pos = s.find_first_of(delims, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos) {
        ss.push_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delims, pos);
        pos = s.find_first_of(delims, lastPos);
        if(ss.size() == (max - 1)) {
            pos = std::string::npos;
        }
    }
}

/***************************/

static inline std::string strip(const std::string& OrgString, const char *Delims)
{

	auto begin_pos = OrgString.find_first_not_of(Delims);
	auto end_pos = OrgString.find_last_not_of(Delims);

	if (begin_pos == std::string::npos || end_pos == std::string::npos) {
		return "";
	}

	return OrgString.substr(begin_pos, end_pos - begin_pos + 1);
}

template <typename Derived, typename Base>
static bool is_class_of(const Base*Node) {
    Derived* DerivedNode = dynamic_cast<Derived *>(const_cast<Base *>(Node));
    return DerivedNode != nullptr;
}

template<typename Iterator>
std::string str_join(Iterator Begin, Iterator End, const std::string& Sep)
{
    std::stringstream ss;

    if (Begin != End) {

        ss << *Begin;
        Begin++;

        while (Begin != End) {
            ss << Sep << *Begin;
            Begin++;
        }
    }

    return ss.str();
}

bool CreatePath(const char *Path);
bool CreatePaths(const char *Path);
std::string SciFloat2Str (float Float);
} // Namespace Utils
#endif /* LCC_UTIL_H */
