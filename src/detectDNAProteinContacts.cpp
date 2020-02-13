#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include"../include/misc/StandardOutput.hpp"
#include<string>
#include<fstream>
#include<iostream>
#include<memory>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
	// check command line arguments
	if (argc != 5) error_output("too much or less arguments");
	cafemol::output_handling::Standard_Output standard_output = cafemol::output_handling::Standard_Output();

	std::array<std::string, 2> input_names = {argv[1], argv[2]};
	std::string output_name = argv[3];
	const float& cutoff = std::stof(argv[4]);

	// default value
	const std::string parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/histone_address";
//	const float cutoff = 10.0;
	const int block_size = 80;
	const int elonged_601_length = 275;

	// suffix
	const std::array<std::string, 2> suffixes = {"dcd", "psf"};

	// check file format
	std::array<std::string, 2> filenames;

	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(input_names, suffixes);
	}

	// -----------------------------------------------------------------------------------------
	// generate instances

	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);


	// -----------------------------------------------------------------------------------------
	// read protein struture file
	standard_output.output_HyphenBlock("", block_size);
	standard_output.output_HyphenBlock("Reading the protein structure information", block_size);

	cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();
	if (chain_data.empty()) error_output("There is no data in this file.");
	std::vector<std::string> atom_name_data = psf_reader->get_AtomNameContainer();
	// read the number of DNA chain residues
	standard_output("Reading the beginning and end residue numbers of each DNA chain.");
	std::vector<std::array<int, 2>> dna_chains_start_ends;

	int total_DNA_atom_num = 0;
	for (const cafemol::psf_data_type::psf_chain_info& chain_info : chain_data) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			dna_chains_start_ends.push_back(chain_start_end);
			total_DNA_atom_num += (chain_start_end[1] - chain_start_end[0] + 1);
		}
	}

	if (dna_chains_start_ends.size() % 2 == 1) error_output("The number of DNA chain is not even.");
	const int dsDNA_number = dna_chains_start_ends.size() / 2;

	std::vector<int> dna_ids;
	for (int dsDNA_idx = 0; dsDNA_idx < dsDNA_number; ++dsDNA_idx) {
		std::array<int, 2> DNA_chainA = dna_chains_start_ends[2 * dsDNA_idx];
		std::array<int, 2> DNA_chainB = dna_chains_start_ends[2 * dsDNA_idx + 1];
		if ((DNA_chainA[1] - DNA_chainA[0]) != (DNA_chainB[1] - DNA_chainB[0])) error_output("This DNA molecules have ssDNA region or The DNA pairing is wrong.");
		for (int dna_id = 0; dna_id < (DNA_chainA[1] - DNA_chainA[0] + 1); ++dna_id) {
			dna_ids.push_back(dna_id + DNA_chainA[0]);
			dna_ids.push_back(DNA_chainB[1] - dna_id);
		}
	}

	standard_output("The number of DNA atoms ->", total_DNA_atom_num, "");

	int unwrapped_num = elonged_601_length - (total_DNA_atom_num + 2) / 6;
	standard_output("The number of processed DNA ->", unwrapped_num);
	standard_output.output_HyphenBlock("", block_size);

	// -----------------------------------------------------------------------------------------
	// read the PDNS residue numbers
	standard_output.output_HyphenBlock("Reading the PDNS protein residues", block_size);
	std::vector<int> pdns_sorted_residues;
	standard_output("Reading the data from ...", parameter_name);
	std::ifstream para_fs(parameter_name, std::ios::in);
	std::string buffer;
	while (std::getline(para_fs, buffer)) {
		pdns_sorted_residues.push_back(total_DNA_atom_num + std::stoi(buffer));
	}
	para_fs.close();
	if (pdns_sorted_residues.empty()) error_output("No data << ", parameter_name, ">>");
	
	standard_output("The number of PDNS contacts ->", pdns_sorted_residues.size());
	standard_output.output_HyphenBlock("", block_size);

	// -----------------------------------------------------------------------------------------
	// find PDNS contact pairs at each MD step
	standard_output.output_HyphenBlock("Finding the PDNS contact pairs at each MD step", block_size);
	std::vector<std::vector<int>> pdns_contact_targets_vec = dcd_analyzer->get_PDNSContactResidues(pdns_sorted_residues, dna_ids, cutoff);
	standard_output.output_HyphenBlock("", block_size);

	// -----------------------------------------------------------------------------------------
	// adjust DNA atom id to DNA base pair
	standard_output.output_HyphenBlock("Adjusting DNA atom ID to DNA bp serial ID", block_size);
	std::vector<std::vector<int>> output_data;

	for (const std::vector<int>& pdns_contact_targets : pdns_contact_targets_vec) {
		output_data.push_back(psf_reader->convert_DNAID2dsDNAResi(pdns_contact_targets));
	}
	standard_output.output_HyphenBlock("", block_size);

	// -----------------------------------------------------------------------------------------
	// making output file
	standard_output.output_HyphenBlock("Making output file", block_size);
	std::ofstream ofs(output_name, std::ios::out);

	for (std::size_t iframe = 0; iframe < output_data.size(); ++iframe) {
		ofs << "Frame " << iframe << std::endl;
		for (std::size_t pdns_idx = 0; pdns_idx < output_data[iframe].size(); ++pdns_idx) {
			if (output_data[iframe][pdns_idx] <= 0) continue;
			ofs << unwrapped_num + output_data[iframe][pdns_idx] << " " << pdns_idx + 1 << std::endl;
		}
	}
	ofs.close();

	standard_output.output_HyphenBlock("", block_size);
	standard_output("Program finished");
	standard_output("");

	return 0;

// old version ---------------------------------------------------------------------------------

//	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
//	// command line arguments
//	if (argc != 5) error_output("too much or less arguments.");
//
//	std::array<std::string, 3> input_names = {argv[1], argv[2], argv[3]};
//	std::string output_name = argv[4];
//
//	// default value
//	const std::string parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/pdns_residue_sort.par";
//	const float cutoff = 10.0;
//	const float angle2rad = acos(-1) / 180.0;
//
//	// suffix
//	const std::array<std::string, 3> suffixes = {"dcd", "psf", "ninfo"};
//	std::array<std::string, 3> filenames;
//
//	// check each file format
//	{
//	    cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();        
//	    filenames = judgement(input_names, suffixes);
//	}
//	
//	// -----------------------------------------------------------------------------------------
//	// calculation
//
//	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
//	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);
//	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(filenames[2]);
//
//
//
//	// read pdns residue file
//	std::cout << "---------Reading the PDNS residues of a nucleosome.-------------------------------" << std::endl;
//
//	std::cout << "Reading the initial pdns contacts" << std::endl;
//
//	cafemol::ninfo_data_type::container_pdns pdns_data = ninfo_reader->read_NinfoPdns();
//	if (pdns_data.empty()) error_output("There is no PDNS contact in this model");
//
//	std::vector<int> pdns_residues;
//	std::vector<std::array<float, 3>> pdns_parameters;
//
//	for (const cafemol::ninfo_data_type::pdns_tuple& pdns_line_data : pdns_data) {
//	    std::array<int, 2> pdns_mps = std::get<2>(pdns_line_data.line_data);
//	    std::array<float, 4> pdns_params = std::get<3>(pdns_line_data.line_data);
//	    pdns_residues.push_back(pdns_mps[0]);
//	    pdns_parameters.push_back({pdns_params[0], angle2rad * pdns_params[1], angle2rad * pdns_params[2]});
//	}
//
//	std::cout << "... Done" << std::endl;
//	std::cout << std::endl;
//	std::cout << "The number of pdns residues -> " << pdns_residues.size() << std::endl;
//	std::cout << "----------------------------------------------------------------------------------" << std::endl;
//
//
//
//	// read protein structure file
//	std::cout << "---------Reading the protein structure information.-------------------------------" << std::endl;
//
//	cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();
//	std::vector<std::string> atom_name_data = psf_reader->get_AtomNameContainer();
//	if (chain_data.empty()) error_output("There is no data in this file");
//
//	// read the number of DNA chain residues
//	std::cout << "Reading the beginning and end residue number of each DNA chain." << std::endl;
//	std::vector<std::array<int, 2>> dna_chains_start_ends;
//
//	for (const cafemol::psf_data_type::psf_chain_info& chain_info : chain_data) {
//	    std::string chain_kind_name = std::get<0>(chain_info);
//	    if (chain_kind_name == "DNA") {
//	        std::array<int, 2> chain_start_end = std::get<1>(chain_info);
//	        std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
//	        dna_chains_start_ends.push_back(chain_start_end);
//	    }
//	}
//	std::cout << "... Done" << std::endl;
//	std::cout << std::endl;
//
//	std::vector<std::array<int, 2>> pdns_phos_sug_ids;
//
//	for (const std::array<int, 2>& dna_chain_start_end : dna_chains_start_ends) {
//	    for (std::size_t idx = dna_chain_start_end[0] - 1; idx <= dna_chain_start_end[1] - 1; ++idx) {
//	        if (atom_name_data[idx] == "DP") pdns_phos_sug_ids.push_back({idx + 1, idx + 2});
//	    }
//	}
//	
//
//	// alignment the residue number.
//	std::cout << "get the residue number of DNA chain." << std::endl;
//	// check whether the number of DNA is even.
//	if (dna_chains_start_ends.size() % 2 == 1) error_output("The number of DNA chain is not even.");
//	const int dsDNA_number = dna_chains_start_ends.size() / 2;
//
//	int DNA_atom_num = 0;
//
//	for (int idx = 0; idx < dsDNA_number; ++idx) {
//	    std::array<int, 2> DNA_chainA = dna_chains_start_ends[2 * idx];
//	    std::array<int, 2> DNA_chainB = dna_chains_start_ends[2 * idx + 1];
//	    if ((DNA_chainA[1] - DNA_chainA[0]) != (DNA_chainB[1] - DNA_chainB[0])) error_output("This DNA molecules have ssDNA region.");
//	    DNA_atom_num += 2 * (DNA_chainA[1] - DNA_chainA[0] + 1);
//	}
//	std::cout << "... Done" << std::endl;
//	std::cout << std::endl;
//	std::cout << "The number of DNA residues -> " << DNA_atom_num << std::endl;
//	std::cout << "----------------------------------------------------------------------------------" << std::endl;
//
//
//
//	// calculate the pdns contact pair at each MD step from DCD file.
//	std::cout << "---------Calculating PDNS pairs at each MD step.----------------------------------" << std::endl;
//
//	std::cout << "reading the file'" << parameter_name << "'" << std::endl;
//	std::ifstream para_ifs(parameter_name, std::ios::in);
//	if (!para_ifs.is_open()) error_output("No such parameter file.");
//	std::string buffer;
//	std::vector<int> pdns_sorted_residues;
//	while (std::getline(para_ifs, buffer)) {
//	    pdns_sorted_residues.push_back(DNA_atom_num + std::stoi(buffer));
//	}
//	para_ifs.close();
//	if (pdns_sorted_residues.empty()) error_output("No data", parameter_name);
//	std::vector<std::vector<std::array<int, 2>>> pdns_contact_pairs = dcd_analyzer->get_PDNSContactResidues(pdns_phos_sug_ids, pdns_residues, pdns_parameters);
//
//	  std::ofstream test("buf2.tes");
//	  for (const std::vector<std::array<int, 2>>& pdns_contact_pair : pdns_contact_pairs) {
//	      for (const std::array<int, 2>& pdns_contact : pdns_contact_pair) {
//	          test << pdns_contact[0] << " " << pdns_contact[1] << std::endl;
//	      }
//	  }
//	  test.close();
//
//	std::cout << "... Done" << std::endl;
//	std::cout << "----------------------------------------------------------------------------------" << std::endl;
//
//	std::vector<std::vector<int>> pdns_contact_dna_ids;
//	std::vector<std::vector<int>> pdns_contact_pro_resis;
//	
//	for (const std::vector<std::array<int, 2>>& pdns_contact_pair : pdns_contact_pairs) {
//	    std::vector<int> dna_id;
//	    std::vector<int> pro_resi;
//	    for (const std::array<int, 2>& pdns_contact : pdns_contact_pair) {
//	        dna_id.push_back(pdns_contact[0]);
//	        int pdns_index = 1;
//	        for (const int& pdns_sorted_residue : pdns_sorted_residues) {
//	            if (pdns_contact[1] == pdns_sorted_residue) break;
//	            ++pdns_index;
//	        }
//	        if (pdns_index >= (pdns_sorted_residues.size() + 1)) error_output("Invalid PDNS residue number exists.");
//	        pro_resi.push_back(pdns_index);
//	    }
//	    pdns_contact_dna_ids.push_back(dna_id);
//	    pdns_contact_pro_resis.push_back(pro_resi);
//	}
//
//	  std::ofstream test("buf.tes");
//	  for (std::size_t idx = 0; idx < pdns_contact_dna_ids.size(); ++idx) {
//	      for (std::size_t jdx = 0; jdx < pdns_contact_dna_ids[idx].size(); ++idx) {
//	          test << pdns_contact_dna_ids[idx][jdx] << " " << pdns_contact_pro_resis[idx][jdx] << std::endl;
//	      }
//	  }
//	  test.close();
//
//	std::vector<std::vector<int>> pdns_contact_dna_bps;
//	for (const std::vector<int>& pdns_contact_dna_id : pdns_contact_dna_ids) {
//	    pdns_contact_dna_bps.push_back(psf_reader->convert_DNAID2dsDNAResi(pdns_contact_dna_id));
//	}
//	if (pdns_contact_pro_resis.size() != pdns_contact_dna_bps.size()) error_output("The number of histone residues is not consistent with that of dna bps");
//
//	// make output files
//	std::cout << "Output PDNS contact pair -> " << output_name << std::endl;
//	std::ofstream ofs(output_name, std::ios::out);
//
//	int iframe = 0;
//	for (std::size_t idx = 0; idx < pdns_contact_dna_bps.size(); ++idx) {
//	    ofs << "Frame " << iframe << std::endl;
//	    if (pdns_contact_pro_resis[idx].size() != pdns_contact_dna_bps[idx].size()) error_output("The number of histone residue is not consistent with that of DNA bp");
//	    for (std::size_t pdns_idx = 0; pdns_idx < pdns_contact_pro_resis[idx].size(); ++pdns_idx) {
//	        ofs << pdns_contact_dna_bps[idx][pdns_idx] << " " << pdns_contact_pro_resis[idx][pdns_idx] << std::endl;
//	    }
//	    ++iframe;
//	}
//
//	std::cout << "PDNS output is over." << std::endl;
//	ofs.close();
//
//	std::cout << "----------------------------------------------------------------------------------" << std::endl;
//	std::cout << "Program finished" << std::endl;
//	std::cout << std::endl;
//
//	return 0;
}
