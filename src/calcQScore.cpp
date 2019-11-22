#include"../include/NINFO/NinfoReader.hpp"
#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/CafemolLibrary.hpp"
#include<iostream>
#include<string>
#include<vector>
#include<array>
#include<memory>
#include<fstream>
#include<cmath>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
	
	if (argc != 5) error_output("too much or less command line arguments.",
								"filename1, filename2, output_name, 'all' or unit-unit");

	// filename
	std::array<std::string, 2> filename_set = {argv[1], argv[2]};
	std::string output_name = argv[3];
	const std::array<std::string, 2> suffix_set = {"dcd", "ninfo"};

	// unit number
	std::string c_units = argv[4];
	std::vector<std::string> c_unit_pairs;
	std::vector<std::array<int, 2>> i_unit_pairs;
	bool is_all_pairs = false;


	if (c_units == "all") {
		is_all_pairs = true;
	}
	else {
		c_unit_pairs = cafemol::library::split_ByHyphen(c_units);

		if (c_unit_pairs[0].size() != c_unit_pairs[1].size()) {
			error_output("The nubmer of left units is not consistent with that of right ones.");
		}
		else if (cafemol::library::isdigit(c_unit_pairs)) {
			i_unit_pairs.resize(c_unit_pairs[0].size());
			for (std::size_t idx = 0; idx < c_unit_pairs[0].size(); ++idx) {
				int i_unit_lhs = static_cast<int>(c_unit_pairs[0][idx] - '0');
				if (i_unit_lhs == 0) {
					i_unit_lhs = 10;
				}
				int i_unit_rhs = static_cast<int>(c_unit_pairs[1][idx] - '0');
				if (i_unit_rhs == 0) {
					i_unit_rhs = 10;
				}
				i_unit_pairs[idx] = {i_unit_lhs, i_unit_rhs};
			}
		}
		else error_output("Invalid input", c_unit_pairs[0], c_unit_pairs[1], "include no digit number.");
	}

	// -----------------------------------------------------------------------------------
	// local constant
	// float cutoff = 6.5; // The cutoff length of native contact in cafemol

	// -----------------------------------------------------------------------------------
	// File open
	std::array<std::string, 2> filenames;
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(filename_set, suffix_set);
	}
	// calculation

	// generate instance
	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(filenames[1]);
	cafemol::ninfo_data_type::container_contact contact_data = ninfo_reader->read_NinfoContact();

	if (contact_data.empty()) error_output("There is no contact in this model.");

	// get imps in unit1 and unit2
	std::cout << "---------Reading the information on native contact pairs--------------------------" << std::endl;
	
	if (!is_all_pairs) {
		std::cout << "Contact Pairs" << std::endl;
		std::cout << " source - target" << std::endl;
		for (const std::array<int, 2>& i_unit_pair : i_unit_pairs) {
			std::cout << i_unit_pair[0] << "	-	" << i_unit_pair[1] << std::endl;
		}
	}
	else {
		std::cout << "Read all contacted pairs" << std::endl;
	}

	std::vector<std::array<int, 2>> contact_pairs;
    std::vector<float> cutoffs;

	if (is_all_pairs) {
		for (const cafemol::ninfo_data_type::contact_tuple& c_tuple : contact_data) {
			std::array<int, 4> pairs = std::get<2>(c_tuple.line_data);
            // std::cout << pairs[0] << " " << pairs[1] << std::endl;
			contact_pairs.push_back({pairs[0], pairs[1]});
            std::array<float, 4> params = std::get<3>(c_tuple.line_data);
            cutoffs.push_back(1.2 * params[0]);
		}
	}
	else {
		for (const cafemol::ninfo_data_type::contact_tuple& c_tuple : contact_data) {
			std::array<int, 2> unit_of_mps = std::get<1>(c_tuple.line_data);
			std::array<int, 4> pairs = std::get<2>(c_tuple.line_data);
            std::array<float, 4> params = std::get<3>(c_tuple.line_data);
			if (cafemol::library::is_containercontains(i_unit_pairs, unit_of_mps)) {
    			contact_pairs.push_back({pairs[0], pairs[1]});
                cutoffs.push_back(1.2 * params[0]);
			}
		}
	}


	
	int ini_contact_num = contact_pairs.size();
	std::cout << std::endl;
	std::cout << "the number of initial native contacts -> " << ini_contact_num << std::endl;

	if (ini_contact_num == 0) {
		std::exit(1);
	}

	std::cout << "----------------------------------------------------------------------------------" << std::endl;
	std::cout << "---------Calculating the number of native contacts in each MD step----------------" << std::endl;

	std::vector<int> contact_numbers = dcd_analyzer->count_NativeContacts(contact_pairs, cutoffs);

	std::cout << "calculate and output qscore to the file '" << output_name << "'" << std::endl;
	std::ofstream ofs(output_name, std::ios::out);
	
	float f_ini_contact_num = static_cast<float>(ini_contact_num);
	for (const int& contact_number : contact_numbers) {
		float f_contact_number = static_cast<float>(contact_number);
		float qscore = f_contact_number / f_ini_contact_num;
		ofs << qscore << std::endl;
        /* for verification */
        // float res = std::round(qscore * 1000) / 1000.00;
		// ofs << res << std::endl;
	}

	std::cout << "----------------------------------------------------------------------------------" << std::endl;
	std::cout << "Program finished" <<  std::endl;
	std::cout << "----------------------------------------------------------------------------------" << std::endl;
	return 0;
}
