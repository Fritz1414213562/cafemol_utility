#include"CSVParser.hpp"
#include<memory>
#include<iostream>

int main() {

	std::unique_ptr<caf::CSVParser> parser = std::make_unique<caf::CSVParser>();

	std::string filename = "charge_H2A.csv";
	parser->load_file(filename);
	parser->write2standard();
	parser->close_file();

	return 0;
}
