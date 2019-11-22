#include"../../include/DCD/DCDAnalyzer.hpp"
#include"../../include/PDB/PDBReader.hpp"
#include"../../include/misc/FileOpenJudge.hpp"
#include<string>
#include<iostream>
#include<fstream>
#include<memory>
#include<array>
#include<vector>


int main(int argc, char *argv[]) {

	// file stream
	std::string dcd_name = argv[1];
	std::string pdb_name = argv[2];
	std::string output_name = argv[3];
	std::string dcd_suffix = "dcd";
	std::string pdb_suffix = "pdb";

	int ARG = 116;
	std::string H3_chain1 = "C";
	std::string H3_chain2 = "G";
	std::string DNAChainA = "A";
	std::string DNAChainB = "B";
	float cutoff = 10.0;

	// file open
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(dcd_name, dcd_suffix);
		judgement->SuffixJudge(pdb_name, pdb_suffix);
	}
	std::unique_ptr<cafemol::DCDAnalyzer> dcd_analyzer = std::make_unique<cafemol::DCDAnalyzer>(dcd_name);
	std::unique_ptr<cafemol::PDBReader> pdb_reader = std::make_unique<cafemol::PDBReader>(pdb_name);
	pdb_reader->read_PDBData();
	
	int dna_atom_size = pdb_reader->get_AtomNumof(DNAChainA);
	dna_atom_size += pdb_reader->get_AtomNumof(DNAChainB);
	int dna_A_size = pdb_reader->get_AtomNumof(DNAChainA);
	std::cout << "dsDNA length -> " << dna_A_size << std::endl;
	int dna_resi_size = pdb_reader->get_ResidueNumof(DNAChainA);

	std::vector<int> dna_ids(dna_atom_size);
	
	dna_ids[0] = 1;
	dna_ids[1] = 2;
	dna_ids[dna_atom_size - 2] = dna_A_size + 2;
	dna_ids[dna_atom_size - 1] = dna_A_size + 1;
	for (int idx = 0; idx < (dna_atom_size - 4) / 6; ++idx) {
		dna_ids[6 * idx + 2] = dna_atom_size - 3 * idx;
		dna_ids[6 * idx + 3] = dna_atom_size - 3 * idx - 1;
		dna_ids[6 * idx + 4] = dna_atom_size - 3 * idx - 2;
		dna_ids[6 * idx + 5] = 3 * idx + 3;
		dna_ids[6 * idx + 6] = 3 * idx + 4;
		dna_ids[6 * idx + 7] = 3 * idx + 5;
	}

	int arg_id1 = pdb_reader->get_AtomIDof(ARG, H3_chain1);
	int arg_id2 = pdb_reader->get_AtomIDof(ARG, H3_chain2);
	std::cout << "H3K116 (C) -> " << arg_id1 << std::endl;
	std::cout << "H3K116 (G) -> " << arg_id2 << std::endl;
//	std::vector<int> ChainC_arg = {arg_id1};
//	std::vector<int> ChainG_arg = {arg_id2};
//
//	std::vector<int> closest_ids1 = dcd_analyzer->getClosestResidue(ChainC_arg, dna_ids, cutoff);
//	std::vector<int> closest_ids2 = dcd_analyzer->getClosestResidue(ChainG_arg, dna_ids, cutoff);
	std::array<int, 2> arginines = {arg_id1, arg_id2};
	std::vector<std::array<int, 2>> closest_ids = dcd_analyzer->getClosestResidue(arginines, dna_ids, cutoff);

	std::ofstream ofs(output_name, std::ios::out);

	std::vector<int> closest_bps;

	for (const std::array<int, 2>& closest_id : closest_ids) {
		int closest_bp;
		int closest_bp1;
		int closest_bp2;

		if (closest_id[0] <= 0) {
			closest_bp1 = -1;
		}
		else if ((closest_id[0] <= dna_A_size) && (0 < closest_id[0])) {
			closest_bp1 = pdb_reader->get_ResidueIDof(closest_id[0]);
		}
		else if (closest_id[0] > dna_A_size) {
			closest_bp1 = pdb_reader->get_ResidueIDof(closest_id[0]);
			closest_bp1 = dna_resi_size - closest_bp1 + 1;
		}

		if (closest_id[1] <= 0) {
			closest_bp2 = -1;
		}
		else if ((closest_id[1] <= dna_A_size) && (0 < closest_id[1])) {
			closest_bp2 = pdb_reader->get_ResidueIDof(closest_id[1]);
		}
		else if (closest_id[1] > dna_A_size) {
			closest_bp2 = pdb_reader->get_ResidueIDof(closest_id[1]);
			closest_bp2 = dna_resi_size - closest_bp2 + 1;
		}

		if ((closest_bp1 < 0) || (closest_bp2 < 0)) {
			closest_bp = -1;
		}
		else {
			closest_bp = (closest_bp1 + closest_bp2) / 2;
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

//	std::vector<int> closest_bps1;
//	std::vector<int> closest_bps2;
//
//	for (const int& closest_id : closest_ids1) {
//		int closest_bp;
//		if (closest_id > dna_A_size) {
//			closest_bp = pdb_reader->get_ResidueIDof(closest_id);
//			closest_bp = dna_resi_size - closest_bp + 1;
//		}
//		else if ((closest_id <= dna_A_size) && (0 < closest_id)) {
//			closest_bp = pdb_reader->get_ResidueIDof(closest_id);
//		}
//		else if (closest_id <= 0) {
//			closest_bp = -1;
//		}
//		closest_bps1.push_back(closest_bp);
//	}
//
//	for (const int& closest_id : closest_ids2) {
//		int closest_bp;
//		if (closest_id > dna_A_size) {
//			closest_bp = pdb_reader->get_ResidueIDof(closest_id);
//			closest_bp = dna_resi_size - closest_bp + 1;
//		}
//		else if ((closest_id <= dna_A_size) && (0 < closest_id)) {
//			closest_bp = pdb_reader->get_ResidueIDof(closest_id);
//		}
//		else if (closest_id <= 0) {
//			closest_bp = -1;
//		}
//		closest_bps2.push_back(closest_bp);
//	}
//
//	for (std::size_t idx = 0; idx < closest_bps1.size(); ++idx) {
//		if ((closest_bps1[idx] < 0) || (closest_bps2[idx] < 0)) {
//			ofs << "nan" << std::endl;
//		}
//		else {
//			int closest_bp1 = closest_bps1[idx];
//			int closest_bp2 = closest_bps2[idx];
//			int closest_bp = (closest_bp1 + closest_bp2) / 2;
//			ofs << closest_bp << std::endl;
//		}
//	}

	ofs.close();

	std::cout << "The program finished" << std::endl;
	std::cout << std::endl;

	return 0;
}
