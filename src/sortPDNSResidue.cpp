#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/NINFO/NinfoReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include<string>
#include<array>
#include<memory>
#include<iostream>
#include<fstream>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();	
	if (argc != 5) error_output("too much or less arguments.");
	// command line arguments
	std::array<std::string, 3> input_names = {argv[1], argv[2], argv[3]};
	std::string output_name = argv[4];
	std::array<std::string, 3> filenames;
	// default value
//	std::string sorted_pdns_output = "pdns_residue_sorted_alongDNA.tes";

	// constant variants

	//// suffix
	std::array<std::string, 3> suffix_names = {"dcd", "psf", "ninfo"};
	//// cutoff length of pdns contact
	const float pdns_cutoff = 8.;

	// check each file format
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(input_names, suffix_names);
	}


//-----------------------------------------------------------------------------------------------
// calculation

	//// generate the instance.

	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(filenames[2]);



	//// read native info file
	std::cout << "-------Reading the native information.--------------------------------------------" << std::endl;
	std::cout << "Reading the initial pdns contact." << std::endl;
	
	cafemol::ninfo_data_type::container_pdns pdns_data = ninfo_reader->read_NinfoPdns();

	std::vector<int> pdns_histone_ids;
	for (const cafemol::ninfo_data_type::pdns_tuple& pdns_line_data : pdns_data) {
		std::array<int, 2> pdns_mps = std::get<2>(pdns_line_data.line_data);
		pdns_histone_ids.push_back(pdns_mps[0]);
	}

	// standard output
	std::cout << "... Done" << std::endl;
	std::cout << std::endl;
	std::cout << "The number of pdns contacts -> " << pdns_histone_ids.size() << std::endl;
	std::cout << "----------------------------------------------------------------------------------" << std::endl;



	//// read protein structure file
	std::cout << "-------Reading the protein structure information.---------------------------------" << std::endl;
	
	cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();

	// read the number of DNA chain residues
	std::cout << "Reading the beginning and end residue number of each DNA chain." << std::endl;
	std::vector<std::array<int, 2>> dna_chains_start_ends;
	std::vector<int> dna_chain_indices;
	int chain_idx = 0;

	for (const cafemol::psf_data_type::psf_chain_info& chain_info : chain_data) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			dna_chains_start_ends.push_back(chain_start_end);
			dna_chain_indices.push_back(chain_idx);
		}
		++chain_idx;
	}
	std::cout << "... Done" << std::endl;
	std::cout << std::endl;

	// alignment the residue number.
	std::cout << "Align the residue number of DNA chain." << std::endl;
	// chech whether the number of DNA is even.
	if (dna_chains_start_ends.size() % 2 == 1) {
		std::cerr << "Error: The number of DNA chain is not even." << std::endl;
		std::exit(1);
	}
	const int dsDNA_number = dna_chains_start_ends.size() / 2;

	std::vector<int> dna_ids;

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
	std::cout << "The number of DNA residue -> " << dna_ids.size() << std::endl;
	std::cout << "----------------------------------------------------------------------------------" << std::endl;



	//// calculate the pdns contacted residue from DCD file.
	std::cout << "-------Sort PDNS contacted residue ID.--------------------------------------------" << std::endl;
	std::cout << "Sort PDNS residue along DNA." << std::endl;
	std::vector<int> sorted_pdns = dcd_analyzer->get_NativeContactsFromIni(dna_ids, pdns_histone_ids, pdns_cutoff);

	std::cout << "----------------------------------------------------------------------------------" << std::endl;

	std::cout << "-------Output Sorted PDNS residues to " << output_name << std::endl;
	std::ofstream ofs(output_name, std::ios::out);
	for (const int& pdns_residue_num : sorted_pdns) {
		ofs << pdns_residue_num << std::endl;
	}
	ofs.close();

	std::cout << "The program finished." << std::endl;
	std::cout << std::endl;

	return 0;
}
