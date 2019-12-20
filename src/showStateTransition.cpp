#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/CafemolLibrary.hpp"
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<array>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != 4) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	std::string contact_file_name = argv[1];
	std::string output_name = argv[2];
	int divide_num = std::stoi(argv[3]);
	std::string suffix = "con";
	std::string fileprefix;

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
		fileprefix = judgement.getFilePrefix();
	}

	// intermediate file
	std::string intermediate_name = fileprefix + "_buffer.im";

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

	para_fs.close();

	// --------------------------------------------------------------------------------------

	if (region_name_list.size() != region_resi_list.size()) eout("Wrong block number");

	int contact_num = 0;
	sout("The number of block_size ->", region_name_list.size());
	for (std::size_t idx = 0; idx < region_name_list.size(); ++idx) {
		contact_num += region_resi_list[idx].size();
		sout("Region: " + region_name_list[idx] + " " + std::to_string(region_resi_list[idx].size()));
	}

	// read contact state file
	
	sout.output_HyphenBlock("Reading the contact state file " + contact_file_name, block_size);
	std::ifstream ifs(contact_file_name, std::ios::in);

	buffer.clear();
	bool is_first = false;
	int frame_read = 0;
//	std::vector<std::vector<std::array<int, 2>>> contacts_pairs;
//	std::vector<std::array<int, 2>> contacts_pair;

	std::vector<char> contact_states(contact_num, 'E');

	sout("Adjust the contact state file to " + intermediate_name);

	std::ofstream intermediate_file(intermediate_name, std::ios::out);
	intermediate_file << contact_num << std::endl;

	while (std::getline(ifs, buffer)) {

		const std::vector<std::string> line_list = cafemol::library::split_String(buffer, delimiter, ignore_char);
		if ((line_list.empty()) || (line_list.size() != 2)) eout("Invalid format in " + contact_file_name);
		if ((line_list[0] == "Frame") && (line_list[1] == "0")) {
			is_first = true;
			continue;
		}
		else if ((line_list[0] == "Frame") && is_first) {
			is_first = false;
			++frame_read;
			continue;
		}
		else if ((line_list[0] == "Frame") && !is_first) {
//			contacts_pairs.push_back(contacts_pair);
			for (const char& state_name : contact_states) {
				intermediate_file << state_name << " ";
			}
			intermediate_file << std::endl;
			++frame_read;
			contact_states.clear();
			contact_states.resize(contact_num, 'E');
			continue;
		}

		if (is_first) continue;

		if (std::stoi(line_list[1]) > contact_num) eout("Histone index is out of range");
		else if (std::stoi(line_list[0]) < divide_num) contact_states[std::stoi(line_list[1]) - 1] = 'W';
		else if (std::stoi(line_list[0]) >= divide_num) contact_states[std::stoi(line_list[1]) - 1] = 'R';
	}
	ifs.close();
	sout(".... Done");

	for (const char& state_name : contact_states) {
		intermediate_file << state_name << " ";
	}
	intermediate_file << std::endl;
	sout("close intermediate file, " + intermediate_name);
	intermediate_file.close();

	sout.output_HyphenBlock("Reading the file, " + intermediate_name, block_size);

	sout("open intermediate file, " + intermediate_name, "");
	std::ifstream intermediate_input(intermediate_name, std::ios::in);
	
	buffer.clear();

	std::vector<std::vector<char>> regions_states;
	// skip the header in the intermediate file
	std::getline(intermediate_input, buffer);

	sout("The contact number in the file, " + intermediate_name + " is  " + buffer);

	sout("Classify the state in each region at each MD steps", "");

	while (std::getline(intermediate_input, buffer)) {
		std::vector<std::string> line_list = cafemol::library::split_String(buffer, delimiter, ignore_char);
		if (line_list.size() != contact_num) eout("The wrong format: The column size is not consistent with contact number, " + std::to_string(contact_num));

		std::vector<char> regions_state(region_resi_list.size());
		for (std::size_t region_idx = 0; region_idx < region_resi_list.size(); ++region_idx) {
			int empty_num = 0;
			int wrapped_num = 0;
			int rewrapping_num = 0;

			for (const int& residue_id : region_resi_list[region_idx]) {
				if (residue_id > contact_num) eout("Histone index is out of range");
				else if (line_list[residue_id - 1] == "E") ++empty_num;
				else if (line_list[residue_id - 1] == "W") ++wrapped_num;
				else if (line_list[residue_id - 1] == "R") ++rewrapping_num;
				else eout("Strange Character exists, '" + line_list[residue_id - 1] + "'");
			}

			if ((empty_num > wrapped_num) && (empty_num > rewrapping_num)) {
				regions_state[region_idx] = 'E';
			}
			else if (wrapped_num > rewrapping_num) regions_state[region_idx] = 'W';
			else regions_state[region_idx] = 'R';
		}
		regions_states.push_back(regions_state);
	}
	intermediate_input.close();

	sout(".... Done");

	sout.output_HyphenBlock("Output the states to " + output_name, block_size);
	std::ofstream ofs(output_name, std::ios::out);

	sout("Output the header", "");
	for (const std::string& region_name : region_name_list) {
		ofs << region_name << " ";
	}
	ofs << std::endl;
	sout(".... Done");

	sout("Output the body", "");
	for (const std::vector<char>& regions_state : regions_states) {
		for (const char& state : regions_state) {
			ofs << state << " ";
		}
		ofs << std::endl;
	}
	sout(".... Done");

	sout("Close the file, " + output_name);
	ofs.close();

//		contacts_pair.push_back({std::stoi(line_list[0]), std::stoi(line_list[1])});

//	if (!contacts_pair.empty()) contacts_pairs.push_back(contacts_pair);
//
//	sout("The frame number -> ", contacts_pairs.size());
//	sout("The last number of frames -> ", frame_read);
//
//	std::vector<std::vector<int>> contacts_states;
//	// contact state
//	// CS = 0 -> empty
//	// CS = 1 -> wrapped
//	// CS = 2 -> rewrapping
//
//	sout("", "calculating the contact states at each steps");
//
//	for (const std::vector<std::array<int, 2>>& pairs_at_step : contacts_pairs) {
//		std::vector<int> each_contact_states;
//		for (const std::vector<int>& region_resi : region_resi_list) {
//
//			int region_size = region_resi.size();
//			std::vector<int> contacts_state;
//			for (const int& resi_idx : region_resi) {
//				for (const std::array<int, 2>& dna_pro_pair : pairs_at_step) {
//					if ((resi_idx == dna_pro_pair[1]) && (dna_pro_pair[0] <= divide_num)) {
//						contacts_state.push_back(0);
//						break;
//					}
//					else if ((resi_idx == dna_pro_pair[1]) && (dna_pro_pair[0] > divide_num)) {
//						contacts_state.push_back(1);
//						break;
//					}
//				}
//			}
//
//			int rewrapping_state_sum = 0;
//			for (const int& contact_state : contacts_state) rewrapping_state_sum += contact_state;
//			int empty_state_num = region_size - contacts_state.size();
//			int wrapped_state_num = contacts_state.size() - rewrapping_state_sum;
//			if ((empty_state_num > rewrapping_state_sum) && (empty_state_num > wrapped_state_num)) {
//				each_contact_states.push_back(0);
//			}
//			else if (rewrapping_state_sum <= wrapped_state_num) each_contact_states.push_back(1);
//			else each_contact_states.push_back(2);
//		}
//		contacts_states.push_back(each_contact_states);
//	}
//	sout(".... Done");
//
//	sout.output_HyphenBlock("Output results to " + output_name, block_size);
//
//	std::ofstream ofs(output_name, std::ios::out);
//
//	for (const std::string& state_name : region_name_list) {
//		ofs << state_name << " ";
//	}
//	ofs << std::endl;
//
//	for (const std::vector<int>& each_contact_states : contacts_states) {
//		for (const int& each_contact_state : each_contact_states) {
//			if (each_contact_state == 0) ofs << "E" << " ";
//			else if (each_contact_state == 1) ofs << "W" << " ";
//			else if (each_contact_state == 2) ofs << "R" << " ";
//			else eout("Invalid value in 'contacts_states'");
//		}
//		ofs << std::endl;
//	}

	sout.output_HyphenBlock("", block_size);
	sout("Program finished", "");

	return 0;
}
