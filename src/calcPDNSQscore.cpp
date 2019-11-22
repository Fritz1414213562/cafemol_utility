#include"../include/NINFO/NinfoReader.hpp"
#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include<iostream>
#include<string>
#include<vector>
#include<array>
#include<memory>
#include<fstream>
#include<cmath>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
	
	if (argc != 6) error_output("too much or less command line arguments.",
								"filename1, filename2, output_name, 'all' or unit-unit");

	// filename
	std::array<std::string, 3> filename_set = {argv[1], argv[2], argv[3]};
	std::string output_name_pdns_proportion = argv[4];
	std::string output_name_pdns_mapping = argv[5];
	const std::array<std::string, 3> suffix_set = {"dcd", "ninfo", "psf"};

	// constant variants
	const float angle2rad = acos(-1) / 180.0;
	const std::string pdns_order_filename = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/pdns_residue_sort.par";

	// File open
	std::array<std::string, 3> filenames;
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(filename_set, suffix_set);
	}
	// calculation

	// generate instance
	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(filenames[1]);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[2]);
	cafemol::ninfo_data_type::container_pdns pdns_data = ninfo_reader->read_NinfoPdns();
	if (pdns_data.empty()) error_output("There is no contact in this model.");

	std::cout << "---------Reading the information on pdns contacts.-------------------------------" << std::endl;
	std::cout << "Reading the initial pdns contact." << std::endl;
	std::vector<int> pdns_histone_ids;
	std::vector<std::array<float, 3>> pdns_parameters;

	for (const cafemol::ninfo_data_type::pdns_tuple& pdns_line_data : pdns_data) {
		std::array<int, 2> pdns_mps = std::get<2>(pdns_line_data.line_data);
		std::array<float, 4> pdns_params = std::get<3>(pdns_line_data.line_data);
		pdns_histone_ids.push_back(pdns_mps[0]);
		pdns_parameters.push_back({pdns_params[0], angle2rad * pdns_params[1], angle2rad * pdns_params[2]});
	}

	std::cout << "... Done" << std::endl;
	std::cout << std::endl;
	std::cout << "The number of pdns contacts -> " << pdns_histone_ids.size() << std::endl;
	std::cout << "---------------------------------------------------------------------------------" << std::endl;

	std::cout << "---------Reading the protein structure information.------------------------------" << std::endl;

	cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();
	std::vector<std::string> atom_name_data = psf_reader->get_AtomNameContainer();
	if (chain_data.empty()) error_output("There is no data in this file");

	
	std::cout << "Reading the beginning and end residue number of each DNA chain." << std::endl;
	std::vector<std::array<int, 2>> dna_chains_start_ends;
	
	for (const cafemol::psf_data_type::psf_chain_info& chain_info : chain_data) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			dna_chains_start_ends.push_back(chain_start_end);
		}
	}
	std::cout << "... Done" << std::endl;
	std::cout << std::endl;

	std::vector<std::array<int, 2>> pdns_phos_sug_ids;

	for (const std::array<int, 2>& dna_chain_start_end : dna_chains_start_ends) {
		for (std::size_t idx = dna_chain_start_end[0] - 1; idx <= dna_chain_start_end[1] - 1; ++idx) {
			if (atom_name_data[idx] == "DP") pdns_phos_sug_ids.push_back({idx + 1, idx + 2});
		}
	}

	std::vector<int> dna_ids;

	const int dsDNA_number = dna_chains_start_ends.size() / 2;
	for (int idx = 0; idx < dsDNA_number; ++idx) {
		std::array<int, 2> DNA_chainA = dna_chains_start_ends[2 * idx];
		std::array<int, 2> DNA_chainB = dna_chains_start_ends[2 * idx + 1];
		if ((DNA_chainA[1] - DNA_chainA[0]) != (DNA_chainB[1] - DNA_chainB[0])) {
			std::cerr << "Error: This DNA molecules have ssDNA region." << std::endl;
			std::exit(1);
		}
		for (int dna_id_idx = 0; dna_id_idx < (DNA_chainA[1] - DNA_chainA[0] + 1); ++dna_id_idx) {
			dna_ids.push_back(dna_id_idx + DNA_chainA[0]);
			dna_ids.push_back(DNA_chainB[1] - dna_id_idx);
		}
	}
	std::cout << "... Done" << std::endl;
	std::cout << std::endl;
	std::cout << "The number of target DNA -> " << pdns_phos_sug_ids.size() << std::endl;
	std::cout << "---------------------------------------------------------------------------------" << std::endl;
	std::cout << "---------Reading DCD trajectory.-------------------------------------------------" << std::endl;
	std::cout << "reading the file '" << pdns_order_filename << "'." << std::endl;
	std::ifstream ifs(pdns_order_filename, std::ios::in);
	int DNA_atom_num = dna_ids.size();
	std::vector<int> pdns_all_ids;
	std::string buffer;
	while (std::getline(ifs, buffer)) {
		pdns_all_ids.push_back(DNA_atom_num + std::stoi(buffer));
	}
	ifs.close();

	std::cout << "getting all PDNS residue number." << std::endl;
	std::cout << "The number -> " << pdns_all_ids.size() << std::endl;
	std::cout << "The first PDNS residue number -> " << pdns_all_ids[0] << std::endl;
	if (pdns_all_ids.empty()) error_output("No data", pdns_order_filename);
	std::vector<std::vector<int>> pdns_contacts = dcd_analyzer->count_PDNSContacts(pdns_phos_sug_ids, pdns_histone_ids, pdns_parameters);
	std::cout << "... Done" << std::endl;
	std::cout << "---------------------------------------------------------------------------------" << std::endl;

	// make output files
	std::cout << "Output PDNS Qscore -> " << output_name_pdns_proportion << std::endl;
	std::ofstream ofs_pdns_qscore(output_name_pdns_proportion, std::ios::out);

	for (const std::vector<int>& pdns_contact : pdns_contacts) {
		float pdns_qscore = static_cast<float>(pdns_contact.size()) / static_cast<float>(pdns_histone_ids.size());
		ofs_pdns_qscore << pdns_qscore << std::endl;
	}
	std::cout << "Qscore calculation is over" << std::endl;
	ofs_pdns_qscore.close();

	std::cout << "Output PDNS Map -> " << output_name_pdns_mapping << std::endl;
	std::ofstream ofs_pdns_map(output_name_pdns_mapping, std::ios::out);

	for (std::size_t iframe = 0; iframe < pdns_contacts.size(); ++iframe) {
		for (std::size_t contact_idx = 0; contact_idx < pdns_all_ids.size(); ++contact_idx) {
			std::size_t is_contacted = 0;
			for (const int& pdns_id : pdns_contacts[iframe]) {
				if (pdns_all_ids[contact_idx] == pdns_id) {
					is_contacted = 1;
					break;
				}
			}
			ofs_pdns_map << iframe << " " << contact_idx << " " << is_contacted << std::endl;
		}
	}
	std::cout << "PDNS Mapping is over." << std::endl;
	ofs_pdns_map.close();

	std::cout << "---------------------------------------------------------------------------------" << std::endl;
	std::cout << "Program finished." << std::endl;
	std::cout << std::endl;

	return 0;
}
