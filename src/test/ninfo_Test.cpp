#include"../include/NINFO/NinfoReader.hpp"
#include<vector>
#include<string>
#include<memory>
#include<iostream>


int main(int argc, char *argv[]) {

	std::string testname = argv[1];
	std::string line_kind_name = "contact";
	
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(testname);
	cafemol::enum_types::ninfo_column contact_num = ninfo_reader->convert_LineKind2Id(line_kind_name);
	std::cout << contact_num << std::endl;
	cafemol::ninfo_data_type::container_contact ninfo_data = ninfo_reader->read_NinfoContact();
	std::cout << ninfo_data.size() << std::endl;

	return 0;
}
