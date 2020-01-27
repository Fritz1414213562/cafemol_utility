#ifndef PARAMETER_FILE_PARSER_HPP
#define PARAMETER_FILE_PARSER_HPP
#include"../misc/UtilFunc.hpp"
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

namespace cafemol {

class Parameter_Parser {

public:

	Parameter_Parser(const std::string& inputfile_name) open_File(inputfile_name);

	~Parameter_Parser() = default;

	void open_File(const std::string& inputfile_name) {
		input_name = inputfile_name;
		open_File();
	}

	void close_File() if (input_file.is_open()) input_file.close();

protected:

	void open_File() {
		close_File();
		input_file.open(input_name);
	}

	std::vector<std::string> split_ParameterLine(const std::string& line) {
		return cafemol::library::split_String(line, delimiter, comment_out_char);
	}

private:

	const char comment_out_char = '*';
	const char delimiter = ' ';
	std::string input_name;
	std::ifstream input_file;

};
}

#endif /* PARAMETER_FILE_PARSER_HPP */
