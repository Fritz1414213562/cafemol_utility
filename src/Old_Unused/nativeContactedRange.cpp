#include"../include/DCD/DCDAnalyzer.hpp"
#include"../include/PDB/PDBReader.hpp"
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
	std::string pdb_name = argv[2];
	std::string output_name = argv[3];
    std::string dcd_suffix = "dcd";
	std::string pdb_suffix = "pdb";
	// ---------------------------------------------------------------------------------
	// constant
	float cutoff = 10.0; // The cutoff length of native contact in cafemol
	std::string His8ChainIDs = "CDEFGHIJ";
	std::string DNAChainID = "AB";
	std::string DNAChainIDA = "A";

	// ---------------------------------------------------------------------------------
	// File open
	{
		std::unique_ptr<cafemol::FileOpenJudge> judgement = std::make_unique<cafemol::FileOpenJudge>();
		judgement->SuffixJudge(dcd_name, dcd_suffix);
		judgement->SuffixJudge(pdb_name, pdb_suffix);
	}

	// ---------------------------------------------------------------------------------
	// calculation

	// generate instance
    std::unique_ptr<cafemol::DCDAnalyzer> dparser =  std::make_unique<cafemol::DCDAnalyzer>(dcd_name);
	std::unique_ptr<cafemol::PDBReader> pparser = std::make_unique<cafemol::PDBReader>(pdb_name);
	pparser->read_PDBData();

	std::cout << "-------Reading the information on DNA from PDB------------------------------------" << std::endl;
	int dna_atom_size = pparser->get_AtomNumof(DNAChainID);
	std::cout << "The number of atoms belonging to dsDNA -> " << dna_atom_size << std::endl;
	int dna_A_size = pparser->get_AtomNumof(DNAChainIDA);
	int dna_resi_size = pparser->get_ResidueNumof(DNAChainIDA);
	std::cout << "The length of base pair -> " << dna_resi_size << std::endl;

	std::vector<int> dna_ids(dna_atom_size);
	// store 2-beads bp
	dna_ids[0] = 1;
	dna_ids[1] = 2;
	dna_ids[dna_atom_size - 2] = dna_A_size + 2;
	dna_ids[dna_atom_size - 1] = dna_A_size + 1;

	if ((dna_atom_size - 4) % 6 != 0) {
		std::cerr << "Error: The number of DNA atoms is wrong" << std::endl;
	}
	for (int idx = 0; idx < (dna_atom_size - 4) / 6; ++idx) {
		dna_ids[6 * idx + 2] = dna_atom_size - 3 * idx;
		dna_ids[6 * idx + 3] = dna_atom_size - 3 * idx - 1;
		dna_ids[6 * idx + 4] = dna_atom_size - 3 * idx - 2;
		dna_ids[6 * idx + 5] = 3 * idx + 3;
		dna_ids[6 * idx + 6] = 3 * idx + 4;
		dna_ids[6 * idx + 7] = 3 * idx + 5;
	}

	std::cout << "-------Reading the information on the Arg in SHL----------------------------------" << std::endl;
	int His8_atom_begin = dna_atom_size + 1;
	int His8_atom_end = pparser->get_AllAtomNum();
	std::cout << "The range of Histone octamer ID -> " << His8_atom_begin << "-" << His8_atom_end << std::endl;

	std::vector<int> His8_ids;
	for (int id = His8_atom_begin; id <= His8_atom_end; ++id) {
		His8_ids.push_back(id);
	}
	
	std::cout << "-------Calculating the 1st base pair contacted to Histone octamer-----------------" << std::endl;
	std::array<std::vector<int>, 2> closest_ids = dparser->get_NativeContactedRange(dna_ids, His8_ids, cutoff);
	std::cout << "-------Caluculation is over-------------------------------------------------------" << std::endl;
	std::cout << "-------Checking the number of atoms in each file. --------------------------------" << std::endl;
	// check atom number written in each file
	int dcd_atom = dparser->get_AtomNumber();
	int pdb_atom = pparser->get_AllAtomNum();
	std::cout << "The number of atoms in PDB file -> " << pdb_atom << std::endl;
	std::cout << "The number of atoms in DCD file -> " << dcd_atom << std::endl;

	if (dcd_atom != pdb_atom) {
		std::cerr << "Error: The number of atoms in PDB is not consistent with that in DCD." << std::endl;
		std::exit(1);
	}

	std::cout << "-------Making an output file------------------------------------------------------" << std::endl;

	std::ofstream ofs(output_name, std::ios::out);

	std::size_t closest_ids_size = closest_ids[0].size();
//	for (std::size_t idx = 0; idx < closest_ids_size; ++idx) {
//		ofs << closest_ids[0][idx] << " " << closest_ids[1][idx] << std::endl;
//	}
	for (std::size_t idx = 0; idx < closest_ids_size; ++idx) {
		int first_closest_id = closest_ids[0][idx];
		int last_closest_id = closest_ids[1][idx];
		int contacted_range = 0;

		if (first_closest_id <= 0) {
//			ofs << contacted_range << std::endl;;
			ofs << "nan nan" << std::endl;
			continue;
		}

		int first_closest_bp = pparser->get_ResidueIDof(first_closest_id);
		int last_closest_bp = pparser->get_ResidueIDof(last_closest_id);

		if ((first_closest_id > dna_A_size) && (last_closest_id > dna_A_size)) {
			first_closest_bp = dna_resi_size - first_closest_bp + 1;
			last_closest_bp = dna_resi_size - last_closest_bp + 1;
		}

		else if ((first_closest_id <= dna_A_size) && (last_closest_id > dna_A_size)) {
			last_closest_bp = dna_resi_size - last_closest_bp + 1;
		}

		else if ((first_closest_id > dna_A_size) && (last_closest_id <= dna_A_size)) {
			first_closest_bp = dna_resi_size - first_closest_bp + 1;
		}

	//	contacted_range = last_closest_bp - first_closest_bp + 1;
	//	ofs << contacted_range << std::endl;
		ofs << first_closest_bp << " " << last_closest_bp  << std::endl;
	}

	std::cout << "The program finished" << std::endl;
	std::cout << std::endl;

    return 0;
}

