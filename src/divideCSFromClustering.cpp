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

	const int& max_argc1 = 7;
	const int& max_argc2 = 8;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if ((argc != max_argc1) && (argc != max_argc2)) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_prefix = argv[1];
	const std::size_t& file_num = std::stoi(argv[2]);
	const std::string& input_suffix = argv[3];
	const std::string& clustering_result = argv[4];
	const std::string& output_prefix = argv[5];
	const std::string& output_suffix = argv[6];
	bool is_set_interval;
	std::size_t interval_length = 1;
	if (argc == 7) is_set_interval = false;
	else if (argc == 8) {
		interval_length = std::stoi(argv[7]);
		is_set_interval = true;
	}

	const std::size_t& BLOCK_SIZE = 90;

	const std::string& s_file_num = argv[2];
	const std::size_t& file_num_digit = s_file_num.size();

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();
	std::unique_ptr<cafemol::ClusteringResultReader> clustering_result_reader = std::make_unique<cafemol::ClusteringResultReader>(clustering_result);
	std::vector<std::vector<std::array<int, 2>>> all_frames;


	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Reading ContactState Data");

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
		if (is_set_interval) {
			contact_state_reader->read_ContactStatesWOIni(all_frames, input_name, interval_length);
		}
		else contact_state_reader->read_ContactStatesWOIni(all_frames, input_name);
	}

	sout("Reading Clustering Result Data");

	std::vector<std::vector<std::size_t>> cluster_ids = clustering_result_reader->read_ClassifiedIDs();

	const std::size_t& cluster_id_digit = std::to_string(cluster_ids.size()).size();

	sout("Divide all frames into clusters");

//	for (std::size_t i_cluster = 0; i_cluster < cluster_ids.size(); ++i_cluster) {
//		const std::size_t& cluster_size_digit = std::to_string(cluster_ids[i_cluster].size()).size();
//		for (std::size_t i_frame = 0; i_frame < cluster_ids[i_cluster].size(); ++i_frame) {
//			const int& frame_index = cluster_ids[i_cluster][i_frame];
//			std::stringstream buffer_filename;
//			buffer_filename << output_prefix
//							<< std::setw(cluster_id_digit)
//							<< std::setfill('0')
//							<< i_cluster + 1
//							<< "-"
//							<< std::setw(cluster_size_digit)
//							<< std::setfill('0')
//							<< i_frame + 1
//							<< output_suffix;
//			std::string output_name(buffer_filename.str());
//			std::ofstream ofs(output_name, std::ios::out);
//			for (const std::array<int, 2>& frame : all_frames[frame_index]) {
//				ofs << frame[0] << " " << frame[1] << std::endl;
//			}
//			ofs.close();
//		}
//	}

	for (std::size_t i_cluster = 0; i_cluster < cluster_ids.size(); ++i_cluster) {
		std::stringstream buffer_filename;
		buffer_filename << output_prefix
						<< std::setw(cluster_id_digit)
						<< std::setfill('0')
						<< i_cluster + 1
						<< output_suffix;
		std::string output_name(buffer_filename.str());
		std::ofstream ofs(output_name, std::ios::out);

		for (std::size_t i_frame = 0; i_frame < cluster_ids[i_cluster].size(); ++i_frame) {
			const std::size_t& frame_index = cluster_ids[i_cluster][i_frame];
			ofs << "Frame " << (i_frame + 1) * interval_length << std::endl;
			for (const std::array<int, 2>& frame : all_frames[frame_index]) {
				ofs << frame[0] << " " << frame[1] << std::endl;
			}
		}
		ofs.close();
	}



	return 0;
}
