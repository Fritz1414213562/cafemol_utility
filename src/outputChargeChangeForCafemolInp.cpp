#include"../include/PAR/ParameterFileReader.hpp"
#include"../include/PDB/PDBReader.hpp"
#include"../include/misc/UtilException.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include<iostream>
#include<vector>
#include<string>
#include<memory>


int main(int argc, char* argv[]) {

	if (argc != 3) {
		std::cerr << "Error: too much or less arguments" << std::endl;
		std::exit(1);
	}

	std::string input_name = argv[1];
	std::string output_name = argv[2];
	const std::string& parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/charge_change.par";

	// open parameter file
	const std::ifstream& ifs(parameter_name, std::ios::in);
	if (!ifs.is_open()) {
		std::cerr << "Error: There is no parameter file. Please check the path" << std::endl;
		std::exit(1);
	}

	// check file prefix
	{
		std::string pdb_suffix = "pdb";
		cafemol::library::FileOpenJudge judgement = cafemol::library::FileOpenJudge();
		judgement(input_name, pdb_suffix);
	}


// ---------------------------------------------------------------------------------------------
	// generate instances
	std::unique_ptr<cafemol::Parameter_Reader> parameter_reader = std::make_unique<cafemol::Parameter_Reader>();
	std::unique_ptr<cafemol::PDBReader> pdb_reader = std::make_unique<cafemol::PDBReader>(input_name);

	// read parameter file
	std::vector<std::string> target_residues;
	std::vector<float> changed_charges;

	std::string line;
	while (std::getline(ifs, line)) {
		std::vector<std::string> line_list;
		try {
			line_list = parameter_reader->split_ParametersOf<2>(line);
		}
		catch (cafemol::library::Util_Exception e) {
			std::cerr << e.what() << std::endl;
			std::exit(1);
		}
		target_residues.push_back(line_list[0]);
		changed_charges.push_back(std::stof(line_list[1]));
	}

	// 


	return 0;
}
