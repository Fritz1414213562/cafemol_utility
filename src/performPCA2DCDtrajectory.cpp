#include<DCD/DCDReader.hpp>
#include<PrincipalComponentAnalysis.hpp>
#include<ErrorMessage.hpp>
#include<FileOpenJudge.hpp>
#include<iostream>
#include<vector>
#include<array>
#include<memory>
#include<string>


int main(int argc, char *argv[]) {

	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	const int& arg_size = 5;
	if (argc != arg_size) eout("too much or less arguments");

	const std::string& input_name = argv[1];
	const std::string& output_name = argv[2];
	const std::string& s_principal_component_size = argv[3];
//	const std::string& s_flag_move2origin = argv[4];
	const std::size_t& principal_component_size = std::stoi(s_principal_component_size);
//	const std::size_t& i_flag_move2origin = std::stoi(s_flag_move2origin);
	const std::string& s_calc_steps = argv[4];
	int&& calc_step_buf = std::stoi(s_calc_steps);
	if (calc_step_buf <= 0) eout("Invalid calculation steps: Negative values");
	const std::size_t& calc_steps = static_cast<std::size_t>(calc_step_buf);

//	if ((i_flag_move2origin != 0) && (i_flag_move2origin != 1)) eout("Invalid flag value, argv[4]");

	// suffix
	const std::string& suffix = "dcd";

	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		judgement(input_name, suffix);
	}

	std::unique_ptr<cafemol::DCDReader> dcd_reader = std::make_unique<cafemol::DCDReader>(input_name);
	cafemol::Trajectory dcd_trajectory = dcd_reader->get_Trajectory();
//	cafemol::EigenTrajectory dcd_trajectory;

	std::unique_ptr<cafemol::analysis::PCA_Performer> PCA = std::make_unique<cafemol::analysis::PCA_Performer>(output_name);
//	std::unique_ptr<cafemol::library::PCA_Performer> PCA = std::make_unique<cafemol::library::PCA_Performer>(output_name);

//	if (i_flag_move2origin == 1) PCA->set_Trajectory2Origin();
//	else {};
	PCA->set_CalculationStep(calc_steps);

	PCA->run(dcd_trajectory, principal_component_size);

	return 0;
}
