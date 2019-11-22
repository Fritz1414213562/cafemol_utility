#include"../../include/NINFO/NinfoReader.hpp"
#include<string>
#include<memory>
#include<iostream>
#include<fstream>

int main(int argc, char *argv[]) {

	std::string ninfo_name = argv[1];
	std::string output_name = "read_pdns_test.txt";
	if (argc > 3) {
		std::cerr << "Error: too much arguments" << std::endl;
		std::exit(1);
	}
	else if (argc == 3) {
		output_name = argv[2];
	}

	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(ninfo_name);
	cafemol::ninfo_data_type::container_pdns pdns_data = ninfo_reader->read_NinfoPdns();

	std::ofstream ofs(output_name, std::ios::out);
	for (const cafemol::ninfo_data_type::pdns_tuple& pdns_line_data : pdns_data) {
		std::array<int, 2> pdns_mps = std::get<2>(pdns_line_data.line_data);
		ofs << pdns_mps[0] << std::endl;
	}

	return 0;
}
