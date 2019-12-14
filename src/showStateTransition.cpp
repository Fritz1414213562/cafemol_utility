#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/CafemolLibrary.hpp"
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<array>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != 5) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	std::string contact_file_name = argv[1];
	std::string output_name = argv[2];
	int divide_num = std::stoi(argv[3]);
	std::string suffix = "con";

	// default value
	//// parameter file
	const std::string parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/histone_devided_region.par";
	const std::string parameter_block_top = "<<<<";
	const std::string parameter_block_bottom = ">>>>";
	//// block size
	const int block_size = 80;
	//// local variants
	const char delimiter = ' ';
	const char ignore_char = '*';
	
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		judgement(contact_file_name, suffix);
	}

	// -----------------------------------------------------------------------------------------
	
	// read paramter file
	sout.output_HyphenBlock("", block_size);
	sout.output_HyphenBlock("Reading the parameter file, " + parameter_name, block_size);

	std::vector<std::string> region_name_list;
	std::vector<std::vector<int>> region_resi_list;
	
	std::ifstream para_fs(parameter_name, std::ios::in);
	std::string buffer;
	std::vector<int> buf_vector;
	bool is_in_block = false;

	while (std::getline(para_fs, buffer)) {
		if (!is_in_block && (buffer.size() < 4)) continue;
		else if (!is_in_block && (buffer.substr(0, 4) == parameter_block_top)) {
			is_in_block = true;
			if (buffer.size() > 6) region_name_list.push_back(buffer.substr(5));
			continue;
		}
		else if (is_in_block && (buffer.size() > 3) && (buffer.substr(0, 4) == parameter_block_bottom)) {
			is_in_block = false;
			region_resi_list.push_back(buf_vector);
			buf_vector.clear();
			continue;
		}
		else if (is_in_block && (cafemol::library::isdigit(buffer))) {
			buf_vector.push_back(std::stoi(buffer));
			continue;
		}
		else eout("Wrong format in " + parameter_name);
	}

	if (region_name_list.size() != region_resi_list.size()) eout("Wrong block number");

	sout("The number of block_size ->", region_name_list.size());
	for (std::size_t idx = 0; idx < region_name_list.size(); ++idx) {
		sout("Region: " + region_name_list[idx] + " " + std::to_string(region_resi_list[idx].size()));
	}

	// read contact state file
	
	sout.output_HyphenBlock("Reading the contact state file" + contact_file_name, block_size);

	std::ifstream ifs(contact_file_name, std::ios::in);
	buffer.clear();
	bool is_first = false;
	int frame_read = 0;

	while (std::getline(ifs, buffer)) {

		const std::vector<std::string> line_list = cafemol::library::split_String(buffer, delimiter, ignore_char);
		if ((line_list.empty()) || (line_list.size() != 2)) eout("Invalid format in " + contact_file_name);
		if ((line_list[0] == "Frame") && (line_list[1] == "0")) {
			is_first = true;
			continue;
		}
		else if (line_list[0] == "Frame") {
			is_first = false;
			++frame_read;
		}
	}

	return 0;
}
