#include<DCD/DCDWriter.hpp>
#include<FileOpenJudge.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<string>
#include<vector>
#include<fstream>
#include<memory>


int main(int argc, char *argv[]) {

	int MAX_ARGUMENTS_NUMBER = 9;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();

	if (argc != MAX_ARGUMENTS_NUMBER) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	std::string input_name = argv[1];
	std::string output_name = argv[2];
	int start_num = std::stoi(argv[3]);
	int save_step = std::stoi(argv[4]);
	int steps = std::stoi(argv[5]);
	int units = std::stoi(argv[6]);
	float delta = std::stof(argv[7]);
	int ver = std::stoi(argv[8]);

	// constant
	std::string dcd_suffix = "dcd";

	// ---------------------------------------------------------------------------------------
	// file open
	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		judgement(input_name, dcd_suffix);
		judgement(output_name, dcd_suffix);
	}

	std::unique_ptr<cafemol::DCDWriter> dcd_writer = std::make_unique<cafemol::DCDWriter>(output_name, input_name);
	dcd_writer->best_fit_convert(start_num, save_step, steps, units, delta, ver);

	return 0;
}
