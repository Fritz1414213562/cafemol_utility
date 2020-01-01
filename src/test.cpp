#include<DCD/DCDWriter.hpp>
#include<ErrorMessage.hpp>
#include<iostream>
#include<string>
#include<memory>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != 3) eout("too much or less arguments");
	std::string inputfilename = argv[1];
	std::string outputfilename = argv[2];

	int i_start = 0;
	int nstep_save = 100;
	int n_tstep = 10000;
	int n_units = 1;
	float delta = 0.3;
	int ver = 24;

	std::unique_ptr<cafemol::DCDWriter> converter = std::make_unique<cafemol::DCDWriter>(outputfilename, inputfilename);

	converter->copy_DCD(i_start, nstep_save, n_tstep, n_units, delta, ver);

	return 0;
}
