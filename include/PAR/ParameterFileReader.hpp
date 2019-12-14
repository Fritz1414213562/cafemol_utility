#ifndef PARAMETER_FILE_READER_HPP
#define PARAMETER_FILE_READER_HPP
#include"ParameterFileParser.hpp"
#include"../util/UtilException.hpp"
#include<string>
#include<vector>

namespace cafemol {

class Parameter_Reader : Parameter_Parser {

public:

	Parameter_Reader(const std::string& inputfile_name) : Parameter_Parser(inputfile_name) {}
	Parameter_Reader() = default;
	~Parameter_Reader() = default;

	template<std::size_t column_size>
	std::vector<std::string> split_ParametersOf(const std::string& line) {
		const std::vector<std::string>& line_list = split_ParameterLine(line);
		if (line_list.size() != column_size) throw library::Util_Exception("The column size is wrong.");
		return line_list;
	}

private:

};
}

#endif /* PARAMETER_FILE_READER_HPP */
