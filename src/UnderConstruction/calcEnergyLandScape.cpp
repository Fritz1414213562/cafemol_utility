#include"../include/PSF/PSFReader.hpp"
#include"../include/misc/FileOpenJudge.hpp"
#include"../include/misc/ErrorMessage.hpp"
#include"../include/misc/StandardOutput.hpp"
#include<string>
#include<fstream>
#include<iostream>
#include<memory>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output error_output = cafemol:error_handling::Error_Output();
	// check command line arguments.
	if (argc != 4) error_output("too much or less arguments");
	cafemol::output_handling::Standard_Output standard_output = cafemol::output_handling::Standard_Output();

	std::array<std::string, 2> input_names = {argv[1], argv[2]};
	std::string output_name = argv[3];

	// default value
	//// parameter file
	const std::string parameter_name = "/home/nagae/cafemol/torus_cafemol/collision_to_torus/utility/para/pdns_residue_sort.par";
	//// block size
	const int block_size = 80;


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

	cafemol::psf_data_

	return 0;
}
