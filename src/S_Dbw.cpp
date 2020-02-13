#include<analysis/S_Dens_bw>
#include<OtherFormat/ClusteringResultReader.hpp>
#include<OtherFormat/ContactStateReader.hpp>
#include<IO/ErrorMessage.hpp>
#include<IO/StandardOutput.hpp>
#include<memory>
#include<vector>
#include<array>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>



int main(int argc, char *argv[]) {

	const int& max_argc = 8;
	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	CafeInLess::IO::Standard_Output sout = CafeInLess::IO::Standard_Output();

	const std::string& input_prefix = argv[1];
	const std::string& s_file_num = argv[2];
	const std::size_t& file_num = std::stoi(s_file_num);
	const std::string& input_suffix = argv[3];
	const std::string& clustering_result = argv[4];
	const std::string& output_name = argv[5];
	const std::size_t& data_frame_x = std::stoi(argv[6]);
	const std::size_t& data_frame_y = std::stoi(argv[7]);

	const std::array<std::size_t, 2>& data_frame_shape = {data_frame_x, data_frame_y};

	const std::size_t& file_num_digit = s_file_num.size();

	const std::size_t& BLOCK_SIZE = 90;

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();
	std::unique_ptr<cafemol::ClusteringResultReader> clustering_result_reader = std::make_unique<cafemol::ClusteringResultReader>(clustering_result);
	std::unique_ptr<CafeInLess::analysis::S_Dens_bw<CafeInLess::analysis::S_DBW_DIST_L2>> s_dens_calculator = std::make_unique<CafeInLess::analysis::S_Dens_bw<CafeInLess::analysis::S_DBW_DIST_L2>>();

	std::vector<std::vector<float>> all_frames;

	sout[BLOCK_SIZE];
	sout("Reading Contact State Data");

	std::ifstream skip_file("SkipNum", std::ios::in);
	std::vector<std::size_t> skip_number;
	if (skip_file.is_open()) {
		sout("The file, SkipNum has read.");
		std::string buffer;
		while (std::getline(skip_file, buffer)) {
			skip_number.push_back(std::stoi(buffer));
		}
	}

	for (std::size_t file_index = 1; file_index <= file_num; ++file_index) {
		bool is_skipped = false;
		for (std::size_t skip_idx = 0; skip_idx < skip_number.size(); ++skip_idx) {
			is_skipped = (is_skipped || (skip_number[skip_idx] == file_index));
		}
		if (is_skipped) {
			sout("File Index " + std::to_string(file_index) + " skipped.");
			continue;
		}
		
		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;
		std::string input_name(buffer_filename.str());
		contact_state_reader->read_ContactStatesWOIni(all_frames, data_frame_shape, input_name);
	}

	sout("Reading Clustering Result Data");
	const std::vector<std::vector<std::size_t>>& cluster_ids = clustering_result_reader->read_ClassifiedIDs();

	sout[BLOCK_SIZE];
	sout("S_Dbw");
	const float& s_dens_result = s_dens_calculator->run(all_frames, cluster_ids);

	sout("Output to " + output_name);
	std::ofstream ofs(output_name, std::ios::out);
	ofs << s_dens_result << std::endl;
	ofs.close();

	return 0;
}
