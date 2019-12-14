#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include"../include/misc/StandardOutput.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
	if (argc != 7) error_output("too much or less arguments");
	cafemol::output_handling::Standard_Output standard_output = cafemol::output_handling::Standard_Output();

	// for file stream
	std::array<std::string, 2> inputfiles = {argv[1], argv[2]};
	std::string output_name = argv[3];
	int i_ChainID = std::stoi(argv[4]);
	int local_ResiID = std::stoi(argv[5]);
	float cutoff = std::stof(argv[6]);
	std::array<std::string, 2> filenames;
	// -------------------------------------------------------------------------------
	// constant
	const std::array<std::string, 2>& suffixes = {"dcd", "pdb"};
	const int block_size = 80;

	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(inputfiles, suffixes);
	}
	// -------------------------------------------------------------------------------
	// generate instances
	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);

	// read protein structure file.
	standard_output.output_HyphenBlock("", block_size);	
	standard_output.output_HyphenBlock("Reading the protein structure information", block_size);

	cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();

	if (chain_date.empty()) error_output("There is no data in this file.");
	standard_output("Reading the beginning and end residue numbers of each DNA chain");
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
	}

	return 0;
}
