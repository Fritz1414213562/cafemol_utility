#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/NINFO/NinfoReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include<string>
#include<memory>
#include<iostream>
#include<fstream>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

	// command line arguments
	std::string dcd_file_name = argv[1];
	std::string psf_file_name = argv[2];
	std::string ninfo_file_name = argv[3];
	std::string output_name = argv[4];
	// default value
//	std::string sorted_pdns_output = "pdns_residue_sorted_alongDNA.tes";

	// constant variants

	//// suffix
	const std::string dcd_suffix = "dcd";
	const std::string psf_suffix = "psf";
	const std::string ninfo_suffix = "ninfo";
	//// cutoff length of pdns contact
	const float cutoff = 10.0;
	const float pdns_cutoff = 12.0;

	// check each file format
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(dcd_file_name, dcd_suffix);
		judgement->SuffixJudge(psf_file_name, psf_suffix);
		judgement->SuffixJudge(ninfo_file_name, ninfo_suffix);
	}


//-----------------------------------------------------------------------------------------------
// calculation

	//// generate the instance.

	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(dcd_file_name);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(psf_file_name);
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(ninfo_file_name);



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
	std::cout << "-------Calculating PDNS contacted residue ID in each snapshot.--------------------" << std::endl;
	std::cout << "Sort PDNS residue along DNA." << std::endl;
	std::vector<int> sorted_pdns = dcd_analyzer->get_NativeContactsFromIni(dna_ids, pdns_histone_ids, pdns_cutoff);

//	std::ofstream test_ofs(sorted_pdns_output, std::ios::out);
//	for (const int& pdns_resi : sorted_pdns) {
//		test_ofs << pdns_resi << std::endl;
//	}
//	test_ofs.close();

	std::cout << "... Done" << std::endl;
	std::cout << "search PDNS contact pairs in each MD step" << std::endl;

	std::vector<std::array<int, 2>> contact_pairs = dcd_analyzer->get_1stNativeContactedResidue(dna_ids, pdns_histone_ids, cutoff);

	std::cout << "... Done" << std::endl;
	std::cout << "----------------------------------------------------------------------------------" << std::endl;
	


	//// make an output file
	std::cout << "-------Output the result to '" << output_name << "'.";
	std::size_t hyphen_len;
	if (output_name.size() > 51) {
		hyphen_len = 0;
	}
	else {
		hyphen_len = 51 - output_name.size();
	}

	for (std::size_t idx = 0; idx < hyphen_len; ++idx) {
		std::cout << "-";
	}
	std::cout << std::endl;

	std::ofstream ofs(output_name, std::ios::out);
//	std::ofstream test("test.test", std::ios::out);


	// calculate the DNA base pair

	std::vector<std::array<int, 2>> residue_num_begin_end = psf_reader->get_ResidueBeginEnds();
	std::vector<std::array<int, 2>> dna_residue_begin_end;
	for (const int& dna_chain_idx : dna_chain_indices) {
		dna_residue_begin_end.push_back(residue_num_begin_end[dna_chain_idx]);
	}

	for (const std::array<int, 2>& contact_pair : contact_pairs) {

		int closest_dna_bp;
		int DNA_chain_num = 0;
		int dsDNA_length = 0;
		bool is_DNA = false;

//		test << contact_pair[0] << std::endl;

		if (contact_pair[0] <= 0) {
			ofs << "nan nan" << std::endl;
			continue;
		}

		// determine which chain DNA atom id belongs to.
		for (const std::array<int, 2>& dna_chain_start_end : dna_chains_start_ends) {
			if ((contact_pair[0] >= dna_chain_start_end[0]) && (contact_pair[0] <= dna_chain_start_end[1])) {
				dsDNA_length = dna_residue_begin_end[DNA_chain_num][1] - dna_residue_begin_end[DNA_chain_num][0] + 1;
				++DNA_chain_num;
				is_DNA = true;
				break;
			}
			++DNA_chain_num;
		}
		if (!is_DNA) {
			std::cerr << "Error: not found the output DNA residue number." << std::endl;
			std::exit(1);
		}

		// adjust the pdns residue
		int pdns_idx = 1;
		for (const int& sorted_pdns_resi : sorted_pdns) {
			if (sorted_pdns_resi == contact_pair[1]) break;
			++pdns_idx;
		}


		if (DNA_chain_num % 2 == 0) {
			closest_dna_bp = psf_reader->get_LocalResidueNumof(contact_pair[0]);
			closest_dna_bp = dsDNA_length - closest_dna_bp + 1;
//			closest_dna_bp = psf_reader->get_ResidueNumof(contact_pair[0]);
		}
		else if (DNA_chain_num % 2 == 1) {
			closest_dna_bp = psf_reader->get_LocalResidueNumof(contact_pair[0]);
		}
		else {
			std::cerr << "Error: Something wrong." << std::endl;
			std::exit(1);
		}

	//	ofs << closest_dna_bp << " " << contact_pair[1] << std::endl;
		ofs << closest_dna_bp << " " << pdns_idx << std::endl;
	}

	ofs.close();
//	test.close();

	std::cout << "----------------------------------------------------------------------------------" << std::endl;

	std::cout << "The program finished." << std::endl;
	std::cout << std::endl;

	return 0;
}
