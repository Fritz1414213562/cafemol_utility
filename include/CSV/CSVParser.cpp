#include"CSVParser.hpp"

void cafemol::CSVParser::load_file(const std::string& filename) {

	close_input_file();
	input_file.open(filename, std::ios::in);
	std::cout << "The file '" << filename << "' has read as csv format."  << std::endl;
}


void cafemol::CSVParser::close_input_file() {

	if (input_file.is_open()) {input_file.close();}
	else {
		std::cout << "The file stream does not have a file." << std::endl;
	}
}


void cafemol::CSVParser::file_read_alert() {

	if (!input_file.is_open()) {
		std::cerr << "Error: The file stream did not read any file." << std::endl;
		std::exit(1);
	}
}


std::vector<std::string> cafemol::CSVParser::split_line(const std::string& line) {

	std::vector<std::string> result;
	std::string buffer = "";

	for (const char& ch : line) {

		std::string buf_ch{ch};

		if ((ch == delimiter) || (ch == end_token)) {
			result.push_back(buffer);
			buffer = "";
		}

		else {
			buffer += buf_ch;
		}
	}

	return result;

}


void cafemol::CSVParser::write2standard() {

	file_read_alert();
	
	std::string line;
	while (std::getline(input_file, line)) {
		line += "\n";

		std::vector<std::string> read_line = split_line(line);
		for (const std::string& buf : read_line) {
			std::cout << buf << " ";
		}
		std::cout << std::endl;
		line = "";
	}


}
