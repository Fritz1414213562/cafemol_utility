#ifndef UTIL_FUNC_HPP
#define UTIL_FUNC_HPP
#include<string>
#include<vector>

namespace cafemol {

namespace library {

std::vector<std::string> split_String(const std::string& line, const char& delim, const char& ignore_char) {

	std::vector<std::string> result;
	std::string buffer;

	for (const char& Char : line) {
		if (((Char == delim) || (Char == ignore_char)) && (buffer.empty())) continue;
		else if ((Char == delim) && (!buffer.empty())) {
			result.push_back(buffer);
			buffer.clear();
		}
		else {
			buffer.push_back(Char);
		}
	}
	if (!buffer.empty()) {
		result.push_back(buffer);
	}
	return result;
}

}
}


#endif /* UTIL_FUNC_HPP */
