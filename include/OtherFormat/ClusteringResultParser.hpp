#ifndef CLUSTERING_RESULT_PARSER_HPP
#define CLUSTERING_RESULT_PARSER_HPP
#include<ErrorMessage.hpp>
#include<UtilFunc.hpp>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

namespace cafemol {

class ClusteringResultParser {

public:

	ClusteringResultParser() = default;
	ClusteringResultParser(const std::string& inputfile_name) : input_name(inputfile_name) {
		open_File();
	}

	~ClusteringResultParser() = default;

	void open_File(const std::string& inputfile_name) {
		input_name = inputfile_name;
		open_File();
	}

	void close_File() {if (input_file.is_open()) input_file.close();}

protected:

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	std::ifstream input_file;

	void open_File() {
		close_File();
		input_file.open(input_name);
		if (!input_file.is_open()) eout("The file, " + input_name + " could not be found.");
	}

	std::vector<std::string> split_Line(const std::string& line) {
		return cafemol::library::split_String(line, delimiter, comment_out_char);
	}

	std::vector<int> split_ClusterIDLine(const std::string& line) {
		return cafemol::library::split_String2Int(line, delimiter);
	}

private:

	const char comment_out_char = '*';
	const char delimiter = ' ';
	std::string input_name;

};
}

#endif /* CLUSTERING_RESULT_PARSER_HPP */
