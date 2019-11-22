#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/NINFO/NinfoReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char *argv[]) {

	// for file stream
    std::string dcd_name = argv[1];
	std::string psf_name = argv[2];
	std::string output_name = argv[3];

	// contant
    std::string dcd_suffix = "dcd";
	std::string psf_suffix = "psf";
	float cutoff = 6.5;

	// ---------------------------------------------------------------------------------
	// File open
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(dcd_name, dcd_suffix);
		judgement->SuffixJudge(psf_name, psf_suffix);
	}

	// ---------------------------------------------------------------------------------
	// calculation

	// generate instance
    std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer =  std::make_unique<cafemol::DCDAnalyzer>(dcd_name);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(psf_name);

	std::cout << "-------------Reading Chain information in PSF file--------------------------------" << std::endl;
	std::cout << "Reading PSF file now." << std::endl;
	std::vector<cafemol::psf_data_type::psf_chain_info> all_chains = psf_reader->get_ChainInfoOfPSF();
	std::cout << "... Done" << std::endl;
	std::cout << "Identifying the Chain name." << std::endl;
	std::vector<std::array<int, 2>> dna_chains_start_ends;
	std::vector<std::array<int, 2>> his8_chains_start_ends;
	
	// store chain ID as the number
	std::vector<int> DNAchainIDs;
	std::vector<int> ProchainIDs;
	int chainID = 1;

	for (const cafemol::psf_data_type::psf_chain_info& chain_info : all_chains) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			dna_chains_start_ends.push_back(chain_start_end);
			DNAchainIDs.push_back(chainID);
		}
		else if (chain_kind_name == "Protein") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			his8_chains_start_ends.push_back(chain_start_end);
			ProchainIDs.push_back(chainID);
		}
		else {
			std::cerr << "Error: Unknown chain name '" << chain_kind_name << "'." << std::endl;
			std::exit(1);
		}
		++chainID;
	}
	std::cout << "... Done" << std::endl;

	std::cout << "----------------------------------------------------------------------------------" << std::endl;

	std::cout << "-------------Calculating chain ID.------------------------------------------------" << std::endl;

	// DNA ids
	std::cout << "DNA ID calculation" << std::endl;
	std::cout << std::endl;

	if (dna_chains_start_ends.size() % 2 == 1) {
		std::cerr << "Error: The number of DNA chains is not even. You don't use dsDNA model." << std::endl;
		std::exit(1);
	}
	int dsDNA_num = dna_chains_start_ends.size() / 2;

	std::vector<int> dna_ids;

	for (int idx = 0; idx < dsDNA_num; ++idx) {
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

	std::cout << "----------------------------------------------------------------------------------" << std::endl;

	std::cout << "-------------Calculating 1st Native Contacted residue ID in Histone chain --------" << std::endl;

	// histone id
	std::vector<int> His8_ids;
	for (const std::array<int, 2>& his8_start_end : his8_chains_start_ends) {
		for (int his8_id = his8_start_end[0]; his8_id <= his8_start_end[1]; ++his8_id) {
			His8_ids.push_back(his8_id);
		}
	}

	std::vector<std::array<int, 2>> contact_pairs = dcd_analyzer->get_1stNativeContactedResidue(dna_ids, His8_ids, cutoff);
	std::vector<int> closest_ids;
	for (const std::array<int, 2>& contact_pair : contact_pairs) {
		int contact_his8_id = contact_pair[1];
		closest_ids.push_back(contact_his8_id);
	}

	std::cout << "-------------Calculation is over.-------------------------------------------------" << std::endl;

	std::cout << "-------------Making an output file -----------------------------------------------" << std::endl;

	std::ofstream ofs(output_name, std::ios::out);

	for (const int& closest_id : closest_ids) {
		if (closest_id <= 0) {
			ofs << "nan" << std::endl;
			continue;
		}

		std::size_t idx = 0;
		for (const std::array<int, 2>& his8_start_end : his8_chains_start_ends) {
			if ((closest_id >= his8_start_end[0]) && (closest_id <= his8_start_end[1])) {
				ofs << ProchainIDs[idx] << std::endl;
				break;
			}
			else ++idx;
		}
		if (idx >= his8_chains_start_ends.size()) {
			std::cerr << "Error: Something wrong! The closest id is out of range." << std::endl;
			std::exit(1);
		}
	}


	std::cout << "Program finished." << std::endl;
	std::cout << std::endl;

    return 0;
}

