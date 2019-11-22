#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char* argv[]) {

	if (argc < 4) {
		std::cerr << "Error: too less arguments" << std::endl;
		std::exit(1);
	}
	else if (argc > 4) {
		std::cerr << "Error: too much arguments" << std::endl;
		std::exit(1);
	}

	// file stream
	std::array<std::string, 2> filename_set = {argv[1], argv[2]};
	std::string output_name = argv[3];
	std::array<std::string, 2> suffix_set = {"dcd", "psf"};
	std::array<std::string, 2> filenames;

	// constant variants
	const float cutoff = 10.0;
	const std::array<float, 3> origin = {0.0, 0.0, 0.0};

	// file open
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(filename_set, suffix_set);
	}

	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(filenames[0]);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);
	
	std::cout << "----------------------------------------------------------------------------------" << std::endl;

	// get topology and chain information
	std::cout << "------Reading the topology of model.---------------" << std::endl;
	cafemol::psf_data_type::container_psf_chain_info all_chains = psf_reader->get_ChainInfoOfPSF();
	std::cout << "... Done" << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "------Getting the Atom IDs at chain begin and end--" << std::endl;

	// store the IDs at beginning and ends.
	std::vector<std::array<int, 2>> dna_chains_start_ends;

	for (const cafemol::psf_data_type::psf_chain_info& chain_info : all_chains) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": ";
			std::cout << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			dna_chains_start_ends.push_back(chain_start_end);
		}
		else if (chain_kind_name == "Protein") continue;
		else {
			std::cerr << "Error: Strange Chain name." << std::endl;
			std::cerr << "		-> " << chain_kind_name << std::endl;
			std::exit(1);
		}
	}
	std::cout << "... Done" << std::endl;

	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "------Aligning the Atom IDs.-----------------------" << std::endl;

	if (dna_chains_start_ends.size() != 2) {
		std::cerr << "Error: too much or less DNA chain." << std::endl;
		std::exit(1);
	}

	std::cout << "DNA ..." << std::endl;

	int dna_A_atom_size = dna_chains_start_ends[0][1] - dna_chains_start_ends[0][0] + 1;
	int dna_B_atom_size = dna_chains_start_ends[1][1] - dna_chains_start_ends[1][0] + 1;

	if (dna_A_atom_size != dna_B_atom_size) {
		std::cerr << "The DNA size of chain A is not consistent with that of chain B." << std::endl;
		std::exit(1);
	}
	int dna_atom_size = dna_A_atom_size + dna_B_atom_size;
	int dna_resi_size = (dna_A_atom_size + 1) / 3;
	std::cout << "The number of DNA atoms -> " << dna_atom_size << std::endl;
	std::cout << "The number of DNA residues -> " << dna_resi_size << std::endl;

	std::vector<int> dna_ids(dna_atom_size);
	
	dna_ids[0] = 1;
	dna_ids[1] = 2;
	dna_ids[dna_atom_size - 2] = dna_A_atom_size + 2;
	dna_ids[dna_atom_size - 1] = dna_A_atom_size + 1;

	for (int idx = 0; idx < (dna_atom_size - 4) / 6; ++idx) {
		dna_ids[6 * idx + 2] = dna_atom_size - 3 * idx;
		dna_ids[6 * idx + 3] = dna_atom_size - 3 * idx - 1;
		dna_ids[6 * idx + 4] = dna_atom_size - 3 * idx - 2;
		dna_ids[6 * idx + 5] = 3 * idx + 3;
		dna_ids[6 * idx + 6] = 3 * idx + 4;
		dna_ids[6 * idx + 7] = 3 * idx + 5;
	}

	std::cout << "---------------------------------------------------" << std::endl;

	std::cout << "------Calculating closest DNA atom ID to H3K116s.--" << std::endl;
	std::vector<int> closest_ids = dcd_analyzer->get_AtomIDsNearPoint(origin, dna_ids, cutoff);
	std::cout << "... Done" << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;

	std::ofstream ofs(output_name, std::ios::out);

	std::vector<int> closest_bps;

	for (const int& closest_id : closest_ids) {
		int closest_bp;

		if (closest_id <= 0) {
			closest_bp = -1;
		}
		else if ((closest_id <= dna_A_atom_size) && (0 < closest_id)) {
			closest_bp = psf_reader->get_LocalResidueNumof(closest_id);
		}
		else if (closest_id > dna_A_atom_size) {
			closest_bp = psf_reader->get_LocalResidueNumof(closest_id);
			closest_bp = dna_resi_size - closest_bp + 1;
		}
		else {
			std::cerr << "Error: Something wrong!" << std::endl;
			std::exit(1);
		}

		closest_bps.push_back(closest_bp);
	}

	for (const int& closest_bp : closest_bps) {
		if (closest_bp <= 0) {
			ofs << "nan" << std::endl;
		}
		else if (closest_bp > 0) {
			ofs << closest_bp << std::endl;
		}
	}

	ofs.close();

	std::cout << "The program finished" << std::endl;
	std::cout << std::endl;
	std::cout << "----------------------------------------------------------------------------------" << std::endl;

	return 0;
}
