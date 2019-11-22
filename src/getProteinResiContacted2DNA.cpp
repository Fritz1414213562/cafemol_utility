#include"../include/NINFO/NinfoReader.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/PDB/PDBReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/CafemolLibrary.hpp"
#include<array>
#include<vector>
#include<fstream>
#include<memory>
#include<iostream>
#include<set>


int main(int argc, char *argv[]) {

	// for file stream
	std::string ninfo_name = argv[1];
	std::string psf_name = argv[2];
	std::string pdb_name = argv[3];
	std::string output_name = argv[4];

	// contact
	std::string ninfo_suffix = "ninfo";
	std::string psf_suffix = "psf";
	std::string pdb_suffix = "pdb";
//	float cutoff = 11.0;
	float cutoff = 10.5;

	// ------------------------------------------------------------------------------
	// File open
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(ninfo_name, ninfo_suffix);
		judgement->SuffixJudge(psf_name, psf_suffix);
		judgement->SuffixJudge(pdb_name, pdb_suffix);
	}

	// ------------------------------------------------------------------------------
	// calculation

	// generate instance
	std::unique_ptr<cafemol::NinfoReader> ninfo_reader = std::make_unique<cafemol::NinfoReader>(ninfo_name);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(psf_name);
	std::unique_ptr<cafemol::PDBReader> pdb_reader = std::make_unique<cafemol::PDBReader>(pdb_name);
	pdb_reader->read_PDBData();


	std::vector<cafemol::psf_data_type::psf_chain_info> all_chains = psf_reader->get_ChainInfoOfPSF();
	cafemol::ninfo_data_type::container_pdns pdns_data = ninfo_reader->read_NinfoPdns();

	std::vector<std::vector<int>> pdns_histone_chains_ids;
	std::vector<int> pdns_histone_ids;
	int current_chain = 0;
	for (const cafemol::ninfo_data_type::pdns_tuple& pdns_line_data : pdns_data) {
		std::array<int, 2> pdns_mp = std::get<2>(pdns_line_data.line_data);
		if (current_chain == 0) {
			current_chain = std::get<1>(pdns_line_data.line_data)[0];
		}
		else if ((std::get<1>(pdns_line_data.line_data)[0] != current_chain)) {
			pdns_histone_chains_ids.push_back(pdns_histone_ids);
			current_chain = std::get<1>(pdns_line_data.line_data)[0];
			pdns_histone_ids.clear();
		}

		pdns_histone_ids.push_back(pdns_mp[1]);
	}
	pdns_histone_chains_ids.push_back(pdns_histone_ids);


	std::vector<std::array<int, 2>> chains_start_ends;
	std::vector<std::vector<std::array<float, 3>>> coordinates_of_chains;
	std::vector<int> dna_chain_indices; // the DNA indices for coordinates_of_chains
	std::vector<int> pro_chain_indices; // the Protein indices for coordinates_of_chains

	int index = 0;
	for (const cafemol::psf_data_type::psf_chain_info& chain_info : all_chains) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			dna_chain_indices.push_back(index);
			chains_start_ends.push_back(chain_start_end);

			std::vector<std::array<float, 3>> coordinates = pdb_reader->get_Coordinateof(chain_start_end[0], chain_start_end[1]);
			coordinates_of_chains.push_back(coordinates);
		}
		else if (chain_kind_name == "Protein") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			pro_chain_indices.push_back(index);
			chains_start_ends.push_back(chain_start_end);

			std::vector<std::array<float, 3>> coordinates = pdb_reader->get_Coordinateof(chain_start_end[0], chain_start_end[1]);
			coordinates_of_chains.push_back(coordinates);
		}
		else {
			std::cerr << "Error: Unknown chain name '" << chain_kind_name << "'." << std::endl;
			std::exit(1);
		}
		++index;
	}

	std::vector<int> closest_ids;

	for (const int& dna_chain_index : dna_chain_indices) {
		std::vector<std::array<float, 3>> dna_coordinates = coordinates_of_chains[dna_chain_index];

		for (const int& pro_chain_index : pro_chain_indices) {
			std::vector<std::array<float, 3>> pro_coordinates = coordinates_of_chains[pro_chain_index];
			std::array<int, 2> pro_chain_start_end = chains_start_ends[pro_chain_index];
			std::vector<int> pdns_indices = pdns_histone_chains_ids[pro_chain_index - 2];
			int pro_chain_start = pro_chain_start_end[0];
//			for (std::size_t idx = 0; idx < pro_coordinates.size(); ++idx) {
//				for (const std::array<float, 3>& dna_vec : dna_coordinates) {
//					float dist_dna_pro = cafemol::library::calc_Distance(pro_coordinates[idx], dna_vec);
//					if (dist_dna_pro <= cutoff) {
//						closest_ids.push_back(pro_chain_start + idx);	
//						break;
//					}
//				}
//			}
			for (std::size_t dna_idx = 0; dna_idx < dna_coordinates.size(); ++dna_idx) {
				int min_id = -1000;
				for (const int& pdns_id : pdns_indices) {
					float dist_dna_pro = cafemol::library::calc_Distance(pro_coordinates[pdns_id - 1], dna_coordinates[dna_idx]);
					if ((dist_dna_pro < cutoff) && !cafemol::library::is_contains<int>(closest_ids, pdns_id + pro_chain_start - 1)) {
						min_id = pdns_id + pro_chain_start - 1;
					}
				}
				if (min_id < 0) continue;
				closest_ids.push_back(min_id);

//				for (std::size_t pro_idx = 0; pro_idx < pro_coordinates.size(); ++pro_idx) {
//					float dist_dna_pro = cafemol::library::calc_Distance(pro_coordinates[pro_idx], dna_coordinates[dna_idx]);
//
//					if (dist_dna_pro < min_dist_dna_pro) {
//						min_dist_dna_pro = dist_dna_pro;
//						min_idx = pro_idx;
//					}
//				}
//				closest_ids.push_back(pro_chain_start + min_idx);
			}
		}
	}


	std::ofstream ofs(output_name, std::ios::out);

	for (const int& closest_id : closest_ids) {
		ofs << closest_id << std::endl;
	}


	return 0;
}
