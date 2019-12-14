#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include"../include/misc/StandardOutput.hpp"
#include"../include/misc/UtilFunc.hpp"
#include<string>
#include<fstream>
#include<iostream>
#include<memory>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
	// check command line arguments.
	if (argc != 4) error_output("too much or less arguments");
	cafemol::output_handling::Standard_Output standard_output = cafemol::output_handling::Standard_Output();

	std::array<std::string, 2> input_names = {argv[1], argv[2]};
	std::string output_name = argv[3];

	// default value
	//// parameter file
	const std::string parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/histone_address.par";
	//// block size
	const int block_size = 80;
	//// for contact state file format
	const char delimiter = ' ';
	const char ignore_char = '*';

	// suffix
	const std::array<std::string, 2> suffixes = {"psf", "con"};

	// check file format
	std::array<std::string, 2> filenames;
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		filenames = judgement(input_names, suffixes);
	}

	// -----------------------------------------------------------------------------------------
	// generate instances

	std::unique_ptr<cafemol::PSFReader> psf_reader = std::make_unique<cafemol::PSFReader>(filenames[0]);

	// read protein structure file
	standard_output.output_HyphenBlock("", block_size);
	standard_output.output_HyphenBlock("Reading the protein structure information", block_size);

	cafemol::psf_data_type::container_psf_chain_info chain_data = psf_reader->get_ChainInfoOfPSF();

	if (chain_data.empty()) error_output("There is no data in this file.");

	// read the number of DNA chain residues
	standard_output("Reading the beginning and end residue numbers of each DNA chain.");

	int total_DNA_atom_num = 0;

	for (const cafemol::psf_data_type::psf_chain_info& chain_info : chain_data) {
		std::string chain_kind_name = std::get<0>(chain_info);
		if (chain_kind_name == "DNA") {
			std::array<int, 2> chain_start_end = std::get<1>(chain_info);
			total_DNA_atom_num += (chain_start_end[1] - chain_start_end[0] + 1);
		}
	}
	total_DNA_atom_num = (total_DNA_atom_num + 2) / 6;

	standard_output("The number of DNA atoms ->", total_DNA_atom_num);
	standard_output.output_HyphenBlock("", block_size);

	// -------------------------------------------------------------------------------------------
	// read the number of PDNS residue
	standard_output.output_HyphenBlock("Reading the PDNS protein residues", block_size);
	
	int total_PDNS_protein_atom_num = 0;
	standard_output("Reading the data from ...", parameter_name);
	std::ifstream para_fs(parameter_name, std::ios::in);

	if (!para_fs.is_open()) error_output("The file, " + parameter_name + ", could not be found.");

	std::string buffer;
	while (std::getline(para_fs, buffer)) {
		++total_PDNS_protein_atom_num;
	}
	para_fs.close();

	standard_output("The number of PDNS contacts ->", total_PDNS_protein_atom_num);
	standard_output.output_HyphenBlock("", block_size);
	
	// -------------------------------------------------------------------------------------------
	// read a contact state file.
	standard_output.output_HyphenBlock("", block_size);
	standard_output("Reading the contact states in ...", filenames[1]);
	
	std::ifstream ifs(filenames[1], std::ios::in);
	if (!ifs.is_open()) error_output("The file, " + filenames[1] + ", could not be found.");

	// result ... DNA x Protein
	std::vector<std::vector<int>> result(total_DNA_atom_num);

	for (std::size_t idx = 0; idx < total_DNA_atom_num; ++idx) {
		result[idx].resize(total_PDNS_protein_atom_num, 0);
	}

	bool is_first = false;
	int frame_read = 0;

	int last_frame_num = 0;

	while (std::getline(ifs, buffer)) {
		const std::vector<std::string> line_list = cafemol::library::split_String(buffer, delimiter, ignore_char);
		if ((line_list.empty()) || (line_list.size() != 2)) error_output("Invalid format in " + filenames[1]);

		if ((line_list[0] == "Frame") && (line_list[1] == "0")) {
			is_first = true;
			continue;
		}
		else if (line_list[0] == "Frame") {
			is_first = false;
			++frame_read;
			last_frame_num = std::stoi(line_list[1]);
			continue;
		}
		
		if (is_first) continue;

		++result[std::stoi(line_list[0]) - 1][std::stoi(line_list[1]) - 1];
	}

	standard_output("The number of frames ->", frame_read);
	standard_output("The last frame number ->", last_frame_num);
	standard_output.output_HyphenBlock("", block_size);

	// making output file
	standard_output.output_HyphenBlock("Making output file", block_size);
	std::ofstream ofs(output_name);
	
	for (const std::vector<int>& pdns_atoms_prob : result) {
		for (const int& pdns_atom_prob : pdns_atoms_prob) {
			float prob = static_cast<float>(pdns_atom_prob) / static_cast<float>(last_frame_num);
			ofs << prob << " ";
		}
		ofs << std::endl;
	}
	ofs.close();

	standard_output.output_HyphenBlock("", block_size);
	standard_output("Program finished");
	standard_output("");

	return 0;
}
