#pragma once
#include<string>
#include<fstream>
#include<iostream>
#include<vector>

namespace cafemol {

class CSVParser {

protected:
	const char delimiter = ',';
	const char end_token = '\n';
	void file_read_alert();
	std::ifstream input_file;
	std::vector<std::string> split_line(const std::string& line);

public:
	CSVParser() = default;
	~CSVParser() = default;
	void load_file(const std::string& filename);
	void close_input_file();
	// for test
	void write2standard();

};

}
