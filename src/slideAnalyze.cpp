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
	int ARG_SHL = 116;
	std::string His8ChainID = "C";
	std::string DNAChainIDA = "A";
	std::string DNAChainIDB = "B";
	float cutoff = 10.0;

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

	std::cout << "-----------------------Reading the information on DNA ----------------------------" << std::endl;
	int dna_atom_size = pparser->get_AtomNumof(DNAChainIDA);
	dna_atom_size += pparser->get_AtomNumof(DNAChainIDB);
	std::cout << "The number of atoms belonging to dsDNA -> " << dna_atom_size << std::endl;
	int dna_A_size = pparser->get_AtomNumof(DNAChainIDA);
	int dna_resi_size = pparser->get_ResidueNumof(DNAChainIDA);
	std::cout << "The length of base pair -> " << dna_resi_size << std::endl;

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

	std::cout << "-----------------------Reading the information on the Arg in SHL -----------------" << std::endl;
	int arg_id = pparser->get_AtomIDof(ARG_SHL, His8ChainID);
	std::cout << "The AtomID of Arg in HistoneH3 -> " << arg_id << std::endl;
	std::vector<int> SHL_ids = {arg_id};
	
	std::cout << "-----------------------Calculating the base pair closest to the Arg in SHL--------" << std::endl;
	std::vector<int> closest_ids = dparser->getClosestResidue(SHL_ids, dna_ids, cutoff);
	std::cout << "-----------------------Caluculation is over---------------------------------------" << std::endl;
	std::cout << "-----------------------Checking the number of atoms in each file -----------------" << std::endl;
	// check atom number written in each file.
	int dcd_atom = dparser->get_AtomNumber();
	int pdb_atom = pparser->get_AllAtomNum();
	std::cout << "The number of atoms in PDB file -> " << pdb_atom << std::endl;
	std::cout << "The number of atoms in DCD file -> " << dcd_atom <<  std::endl;

	if (dcd_atom != pdb_atom) {
		std::cerr << "Error: The number of atoms in PDB file is not consistent with that in DCD." << std::endl;
		std::exit(1);
	}

	std::cout << "-----------------------Making an output file--------------------------------------" << std::endl;

	std::ofstream ofs(output_name, std::ios::out);

	for (const int& closest_id : closest_ids) {
		int closest_bp;
		if (closest_id > dna_A_size) {
			closest_bp = pparser->get_ResidueIDof(closest_id);
			closest_bp = dna_resi_size - closest_bp + 1;
			ofs << closest_bp << std::endl;
		}
		else if ((closest_id <= dna_A_size) && (0 < closest_id)) {
			closest_bp = pparser->get_ResidueIDof(closest_id);
			ofs << closest_bp << std::endl;
		}

		else if (closest_id <= 0) {
			ofs << "nan" << std::endl;
		}
	}

	std::cout << "The program finished" << std::endl;
	std::cout << std::endl;

    return 0;
}

