#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<OtherFormat/ContactStateReader.hpp>
#include<fstream>
#include<string>
#include<memory>
#include<vector>
#include<array>


int main(int argc, char *argv[]) {

	const int& max_argc = 4;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_name = argv[1];
	const std::string& output_name = argv[2];
	const float& cutoff_len = std::stof(argv[3]);

	const std::size_t& BLOCK_SIZE = 90;

	const std::string& suffix = "dcd";

	{
		cafemol::FileOpenJudge judgement = cafemol::FileOpenJudge();
		judgement(input_name, suffix);
	}

	// ------------------------------------------------------------------------------------------
	// calculation
	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Calculation starts.", "");

	// trajectory
	std::unique_ptr<cafemol::DCDReader> trajectory_reader = std::make_unique<cafemol::DCDReader>(input_name);
	//  contact search
	std::unique_ptr<cafemol::ContactSearcher> contact_search = std::make_unique<cafemol::ContactSearcher>(cutoff_len);

	// reading a trajectory
	sout("Reading a trajectory");
	const cafemol::Trajectory& dcd_trajectory = trajectory_reader->get_Trajectory();
	sout(".... Done", "");

	// contact search
	sout("Searching contacted pairs, and make a contact map");
	if ((dcd_trajectory[0][0].size() != dcd_trajectory[0][1].size()) || (dcd_trajectory[0][0].size() != dcd_trajectory[0][2].size())) eout("The degree of freedom is not consistent among xyz.");
	const std::size_t& atom_size = dcd_trajectory[0][0].size();
	const float frame_size = dcd_trajectory.size();

	std::vector<std::vector<float>> result(atom_size);
	for (std::size_t idx = 0; idx < atom_size; ++idx) {
		result[idx].resize(atom_size, 0);
	}
	
	for (const std::array<std::vector<float>, 3>& snapshot : dcd_trajectory) {
		const std::vector<std::vector<float>>& distance_map = contact_search->get_DistanceMatrix(snapshot);
		for (std::size_t i_row = 0; i_row < atom_size - 1; ++i_row) {
			for (std::size_t i_col = i_row; i_col < atom_size; ++i_col) {
				if (distance_map[i_row][i_col] < 0) {}
				else if (i_row != i_col) {
					result[i_row][i_col] += 1;
					result[i_col][i_row] += 1;
				}
				else {
					result[i_row][i_col] += 1;
				}
			}
		}
	}

	sout(".... Done");
	sout("Calculation ends.", "");

	sout.output_HyphenBlock("Output the calculation results to " + output_name, BLOCK_SIZE);
	std::ofstream ofs(output_name, std::ios::out);

	for (std::size_t i_row = 0; i_row < atom_size; ++i_row) {
		for (std::size_t i_col = 0; i_col < atom_size; ++i_col) {
			result[i_row][i_col] /= frame_size;
			ofs << i_row << " " << i_col << " " << result[i_row][i_col] << std::endl;
		}
	}
	ofs.close();

	sout("Program finished", "");


	return 0;
}
