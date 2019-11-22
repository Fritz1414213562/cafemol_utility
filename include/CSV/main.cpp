#include"makeRESPACData.hpp"
#include<memory>
#include<iostream>
#include<string>

int main(int argc, char* argv[]) {

	if (argc < 4) {
		std::cerr << "Error: The number of arguments is 1." << std::endl;
		std::exit(1);
	}

	// substitute the argument.
	std::string output_name = argv[argc - 1];
	std::string line_prefix = argv[1];

	// file process
	std::unique_ptr<caf::DataMaker> parser = std::make_unique<caf::DataMaker>();
	parser->setLinePrefix(line_prefix);
	for (std::size_t i = 2; i < argc - 1; ++i) {
		std::string filename = argv[i];
		// read file once
		parser->load_file(filename);
		parser->makeChargeData<int>(output_name, 1, 3);
		parser->close_input_file();
		// read file twice
		parser->load_file(filename);
		parser->makeChargeData<int>(output_name, 5, 7);
		parser->close_input_file();
	}

	// close output file
	parser->close_output_file();

	return 0;
}
