#include<OtherFormat/ContactStateReader.hpp>
#include<OtherFormat/ClusteringResultReader.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<sstream>
#include<iomanip>
#include<fstream>
#include<memory>
#include<array>
#include<vector>
#include<string>


int main(int argc, char *argv[]) {

	const int& max_argc = 7;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_prefix = argv[1];
	const std::size_t& file_num = std::stoi(argv[2]);
	const std::string& input_suffix = argv[3];
	const std::string& clustering_result = argv[4];
	const std::string& output_prefix = argv[5];
	const std::string& output_suffix = argv[6];

	const std::size_t& BLOCK_SIZE = 90;

	const std::string& s_file_num = argv[2];
	const std::size_t& file_num_digit = s_file_num.size();

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();
	std::unique_ptr<cafemol::ClusteringResultReader> clustering_result_reader = std::make_unique<cafemol::ClusteringResultReader>(clustering_result);
	std::vector<std::vector<std::array<int, 2>>> all_frames;

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Reading ContactState Data");

	for (std::size_t file_index = 1; file_index <= file_num; ++file_index) {
		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;
		std::string input_name(buffer_filename.str());
		contact_state_reader->read_ContactStatesWOIni(all_frames, input_name);
	}

	sout("Reading Clustering Result Data");

	std::vector<std::vector<int>> cluster_ids = clustering_result_reader->read_ClassifiedIDs();

	const std::size_t& cluster_id_digit = std::to_string(cluster_ids.size()).size();

	sout("Divide all frames into clusters");

	for (std::size_t i_cluster = 0; i_cluster < cluster_ids.size(); ++i_cluster) {
		const std::size_t& cluster_size_digit = std::to_string(cluster_ids[i_cluster].size()).size();
		for (std::size_t i_frame = 0; i_frame < cluster_ids[i_cluster].size(); ++i_frame) {
			const int& frame_index = cluster_ids[i_cluster][i_frame];
			std::stringstream buffer_filename;
			buffer_filename << output_prefix
							<< std::setw(cluster_id_digit)
							<< std::setfill('0')
							<< i_cluster + 1
							<< "-"
							<< std::setw(cluster_size_digit)
							<< std::setfill('0')
							<< i_frame + 1
							<< output_suffix;
			std::string output_name(buffer_filename.str());
			std::ofstream ofs(output_name, std::ios::out);
			for (const std::array<int, 2>& frame : all_frames[frame_index]) {
				ofs << frame[0] << " " << frame[1] << std::endl;
			}
			ofs.close();
		}
	}


	return 0;
}
