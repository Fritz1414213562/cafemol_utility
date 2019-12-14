#include"../include/PDB/PDBReader.hpp"
#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/CafemolLibrary.hpp"
#include<array>
#include<vector>
#include<fstream>
#include<memory>
#include<iostream>
#include<set>


int main(int argc, char *argv[]) {
	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
	if (argc != 4) error_output("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	// for file stream
	std::array<std::string, 2> input_names = {argv[1], argv[2]};
	std::string output_name = argv[3];
	// paramter file
	std::string paramter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/sort_pdns_residue.par";

	std::array<std::string, 2> filenames;

	// constant
	const float cutoff_ele = 5.5;
	const float cutoff_ele_sqr = cutoff_ele * cutoff_ele;
	const float cutoff_pdns = 10.5;
	const float cutoff_pdns_sqr = cutoff_pdns * cutoff_pdns;
	const int block_size = 80;
	const std::array<std::string, 2> suffixes = {"pdb", "psf"};

	// ------------------------------------------------------------------------------------------
	// file open
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(input_names, suffixes);
	}

	// ------------------------------------------------------------------------------------------
	// calculation

	// generate instances
	std::unique_ptr<cafemol::PDBReader> pdb_reader = std::make_unique<cafemol::PDBReader>(filenames[0]);
	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[1]);

	// read PDB
	sout.output_HyphenBlock("", block_size);

	sout("Reading the PDB....");
	pdb_reader->read_PDBData();

	// read a protein structure file
	sout("Reading the protein structure file....");
	cafemol::psf_data_type::container_psf_chain_info all_chains = psf_reader->get_ChainInfoOfPSF();
	if (all_chains.empty()) error_output("There is no data in the file, " + filenames[1] + ".");
	std::vector<std::string> atom_name_data = psf_reader->get_AtomNameContainer();
	std::vector<std::string> resi_name_data = psf_reader->get_ResidueNameContainer();

	std::vector<int> phosphate_ids_unsort;
	std::vector<int> ele_ids;
	int total_DNA_atom_num = 0;

	for (const cafemol::psf_data_type::psf_chain_info chain_info : all_chains) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			int DNA_atom_num = (chain_start_end[1] - chain_start_end[0] + 1);
			total_DNA_atom_num += DNA_atom_num;
			for (std::size_t idx = chain_start_end[0] - 1; idx < chain_start_end[1]; ++idx) {
				if (atom_name_data[idx] == "DP") phosphate_ids_unsort.push_back(idx + 1);
			//	phosphate_ids_unsort.push_back(idx + 1);
			}
		}
		else if (chain_kind_name == "Protein") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			std::cout << chain_kind_name << ": " << chain_start_end[0] << "-" << chain_start_end[1] << std::endl;
			for (std::size_t idx = chain_start_end[0] - 1; idx < chain_start_end[1]; ++idx) {
				if (resi_name_data[idx] == "LYS") ele_ids.push_back(idx + 1);
				else if (resi_name_data[idx] == "ARG") ele_ids.push_back(idx + 1);
				else continue;
			}
		}
		else {
			std::cout << "Warning: Strange chain name '" << chain_kind_name << "'" << std::endl;
			continue;
		}
	}

	std::vector<int> phosphate_ids;
	int ds_DNA_phos_size = phosphate_ids_unsort.size() / 2;
	for (std::size_t idx = 0; idx < ds_DNA_phos_size; ++idx) {
		phosphate_ids.push_back(phosphate_ids_unsort[idx]);
		phosphate_ids.push_back(phosphate_ids_unsort[phosphate_ids_unsort.size() - 1 - idx]);
	}

//	for (const int& phos_id : phosphate_ids) {
//		std::cout << phos_id << std::endl;
//	}


	sout("... Done");
	sout.output_HyphenBlock("", block_size);

	// read paramter file.
	sout.output_HyphenBlock("Reading pdns residues from " + paramter_name, block_size);
	std::ifstream ifs(paramter_name, std::ios::in);
	std::string buffer;
	std::vector<int> pdns_ids;
	
	while (std::getline(ifs, buffer)) {
		pdns_ids.push_back(std::stoi(buffer) + total_DNA_atom_num);
	}
	ifs.close();
	sout.output_HyphenBlock("", block_size);

	// read coordinates and detect contacts

	sout.output_HyphenBlock("Reading Coordinates from PDB", block_size);
	sout("calculating the distance between DNA and pdns or cation residue");
	std::vector<int> sorted_protein_ids;
	std::vector<int> pdns_log;
	std::vector<int> ele_log;
	int pdns_detected = 0;
	int ele_detected = 0;

//	for (const int& phosphate_id : phosphate_ids) {
//		const std::array<float, 3> phosphate_vector = pdb_reader->get_Coordinateof(phosphate_id);
//		float min_pdns_dist_sqr = cutoff_pdns_sqr;
//		float min_ele_dist_sqr = cutoff_ele_sqr;
//		int closest_pdns_id = 0;
//		int closest_ele_id = 0;
//
//		for (const int& pdns_id : pdns_ids) {
//			float pdns_dist_sqr = 0.0;
//			const std::array<float, 3> pdns_vector = pdb_reader->get_Coordinateof(pdns_id); 
//			for (std::size_t idim = 0; idim < 3; ++idim) {
//				pdns_dist_sqr += (pdns_vector[idim] - phosphate_vector[idim]) *
//								 (pdns_vector[idim] - phosphate_vector[idim]);
//			}
//			if (pdns_dist_sqr > min_pdns_dist_sqr) continue;
//			min_pdns_dist_sqr = pdns_dist_sqr;
//			closest_pdns_id = pdns_id;
//		}
//		if ((closest_pdns_id != 0) && (!cafemol::library::is_contains<int>(pdns_log, closest_pdns_id))) {
//			pdns_log.push_back(closest_pdns_id);
//		}
//		pdns_detected = pdns_log.size();
//
//		for (const int& ele_id : ele_ids) {
//			float ele_dist_sqr = 0.0;
//			const std::array<float, 3> ele_vector = pdb_reader->get_Coordinateof(ele_id);
//			for (std::size_t idim = 0; idim < 3; ++idim) {
//				ele_dist_sqr += (ele_vector[idim] - phosphate_vector[idim]) *
//								(ele_vector[idim] - phosphate_vector[idim]);
//			}
//			if (ele_dist_sqr > min_ele_dist_sqr) continue;
//			min_ele_dist_sqr = ele_dist_sqr;
//			closest_ele_id = ele_id;
//		}
//		if ((closest_ele_id != 0) && (!cafemol::library::is_contains<int>(ele_log, closest_ele_id))) {
//			ele_log.push_back(closest_ele_id);
//		}
//		ele_detected = ele_log.size();
//
//		if ((closest_pdns_id == 0) || (cafemol::library::is_contains<int>(sorted_protein_ids, closest_pdns_id))) {}
//		else {
//			sorted_protein_ids.push_back(closest_pdns_id);
//		}
//
//		if ((closest_ele_id == 0) || (cafemol::library::is_contains<int>(sorted_protein_ids, closest_ele_id))) {}
//		else {
//			sorted_protein_ids.push_back(closest_ele_id);
//		}
//	
//	}
//	sout("The number of detected PDNS ->", pdns_detected);
//	sout("The number of detected ELE ->", ele_detected);
//	sout("... Done");
//	sout("The number of interacted residues ->", sorted_protein_ids.size());
//
//	sout.output_HyphenBlock("output results to " + output_name, block_size);
//
//	std::ofstream ofs(output_name, std::ios::out);
//	for (const int& sorted_protein_id : sorted_protein_ids) {
//		ofs << sorted_protein_id << std::endl;
//	}
//	sout("... Done");
//	sout("close the file, " + output_name);
//	ofs.close();
//
//	sout("", "Program finished");

//	std::vector<std::array<int, 2>> dna_pro_pairs;
	std::vector<int> dna_sort;
	std::vector<int> pro_sort;
	sout("calculating the relative distance between PDNS residue and phosphate");
	for (const int& pdns_id : pdns_ids) {
		const std::array<float, 3> pdns_vector = pdb_reader->get_Coordinateof(pdns_id);
		std::array<int, 2> dna_pro_pair = {0, pdns_id};
		int pdns_min_dist_sqr = cutoff_pdns_sqr;
		int phos_idx = 1;
		for (const int& phosphate_id : phosphate_ids) {
			float pdns_dist_sqr = 0.0;
			const std::array<float, 3> phosphate_vector = pdb_reader->get_Coordinateof(phosphate_id);
			for (std::size_t idim = 0; idim < 3; ++idim) {
				pdns_dist_sqr += (phosphate_vector[idim] - pdns_vector[idim]) *
								 (phosphate_vector[idim] - pdns_vector[idim]);
			}
			if (pdns_dist_sqr > pdns_min_dist_sqr) {
				++phos_idx;
				continue;
			}
			pdns_min_dist_sqr = pdns_dist_sqr;
			dna_pro_pair[0] = phos_idx;
			++phos_idx;
		}
		if (dna_pro_pair[0] == 0) continue;
	//	dna_pro_pairs.push_back(dna_pro_pair);
		dna_sort.push_back(dna_pro_pair[0]);
		pro_sort.push_back(dna_pro_pair[1]);
	}
	
	sout("... Done");
//	int detected_pdns_num = dna_pro_pairs.size();
	int detected_pdns_num = pro_sort.size();
	sout("The number of detected PDNS residues -> ", detected_pdns_num);

	sout("calculating the relative distance between cation residue and phosphate");
	for (const int& ele_id : ele_ids) {
		const std::array<float, 3> ele_vector = pdb_reader->get_Coordinateof(ele_id);
		std::array<int, 2> dna_pro_pair = {0, ele_id};
		int ele_min_dist_sqr = cutoff_ele_sqr;
		int phos_idx = 1;
		for (const int& phosphate_id : phosphate_ids) {
			float ele_dist_sqr = 0.0;
			const std::array<float, 3> phosphate_vector = pdb_reader->get_Coordinateof(phosphate_id);
			for (std::size_t idim = 0; idim < 3; ++idim) {
				ele_dist_sqr += (phosphate_vector[idim] - ele_vector[idim]) *
								(phosphate_vector[idim] - ele_vector[idim]);
			}
			if (ele_dist_sqr > ele_min_dist_sqr) {
				++phos_idx;
				continue;
			}
			ele_min_dist_sqr = ele_dist_sqr;
			dna_pro_pair[0] = phos_idx;
			++phos_idx;
		}
		if (dna_pro_pair[0] == 0) continue;
		else if (cafemol::library::is_contains<int>(pro_sort, dna_pro_pair[1])) continue;
	//	dna_pro_pairs.push_back(dna_pro_pair);
		dna_sort.push_back(dna_pro_pair[0]);
		pro_sort.push_back(dna_pro_pair[1]);
	}
	sout("... Done");
	if (dna_sort.size() != pro_sort.size()) error_output("DNA and Protein sort was wrong.");
	sout("The number of detected electrostatic interacted residues -> ", dna_sort.size() - detected_pdns_num);
	sout("The number of total interacted residues -> ", dna_sort.size());

	// sort protein residue;
	sout.output_HyphenBlock("", block_size);
	sout("sort the residues along phosphate ids order");
	std::vector<int> result;
	for (std::size_t idx = 1; idx <= phosphate_ids.size(); ++idx) {
		for (std::size_t pair_idx = 0; pair_idx < dna_sort.size(); ++pair_idx) {
			if (idx == dna_sort[pair_idx]) {
				result.push_back(pro_sort[pair_idx]);
			}
		}
	}
	sout("... Done");

	sout.output_HyphenBlock("output results to " + output_name, block_size);
  
	std::ofstream ofs(output_name, std::ios::out);
	for (const int& res : result) {
		ofs << res << std::endl;
	}
	sout("... Done");
	sout("close the file, " + output_name);
	ofs.close();
  
	sout("", "Program finished");

  return 0;
}

