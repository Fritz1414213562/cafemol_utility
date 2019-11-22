#include"CSVParser.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<cstdarg>

namespace cafemol {

class DataMaker : public CSVParser {

private:
	std::ofstream output_file;
	std::string line_prefix;

public:
	DataMaker() = default;
	~DataMaker() = default;

	void setLinePrefix(const std::string& prefix_name) {line_prefix = prefix_name;}

	template<class... C>
	void makeChargeData(const std::string& filename, C... columns) {

		file_read_alert();
		open_output_file(filename);

		std::string line;
		while (std::getline(input_file, line)) {
			line += "\n";

			std::vector<std::string> read_line = split_line(line);

			output_file << line_prefix << " ";
			for (const int& column : std::initializer_list<int>{ columns... }) {
				output_file << read_line[column] << " ";
			}
			output_file << std::endl;
			line = "";
		}
	}
	
	void open_output_file(const std::string& filename) {
		if (!output_file.is_open()) output_file.open(filename, std::ios::out);
	}

	void close_output_file() {if (output_file.is_open()) output_file.close();}

};

}
