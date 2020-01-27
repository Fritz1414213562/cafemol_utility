#include<OtherFormat/ContactStateReader.hpp>
#include<ContactStateClustering.hpp>
#include<StandardOutput.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<string>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<memory>
#include<array>
#include<vector>
#include<chrono>


int main(int argc, char *argv[]) {

	const int& max_argc = 10;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	// file stream
	const std::string input_prefix = argv[1];
	const std::size_t file_num = std::stoi(argv[2]);
	const std::string input_suffix = argv[3];
	const std::string output_name = argv[4];

	// analysis set-up
	const std::size_t max_cluster_num = std::stoi(argv[5]);
	const std::size_t max_iteration_num = std::stoi(argv[6]);
	const float cutoff = std::stof(argv[7]);
	const std::size_t n_seed = std::stoi(argv[8]);
	const float bin_w = std::stoi(argv[9]);

	const std::size_t BLOCK_SIZE = 90;

	// check a digit of file number
	const std::string s_file_num = argv[2];
	const std::size_t& file_num_digit = s_file_num.size();


	// check Skip Number
	const std::string& skip_filename = "SkipNum";
	std::ifstream skip_file(skip_filename, std::ios::in);
	std::vector<std::size_t> skip_number;
	if (skip_file.is_open()) {
		std::string buffer;
		while (std::getline(skip_file, buffer)) {
			skip_number.push_back(std::stoi(buffer));
		}
	}
	skip_file.close();


// body

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();
	std::unique_ptr<cafemol::analysis::ContactStateClustering> analyzer = std::make_unique<cafemol::analysis::ContactStateClustering>(max_cluster_num, max_iteration_num, cutoff, n_seed, bin_w);

	std::vector<std::vector<std::array<float, 3>>> all_frames;

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Reading Contact State Data");
	

	for (std::size_t file_index = 1; file_index <= file_num; ++file_index) {
		bool is_skipped = false;
		for (std::size_t skip_idx = 0; skip_idx < skip_number.size(); ++skip_idx) {
			is_skipped = (is_skipped || (skip_number[skip_idx] == file_index));
		}
		if (is_skipped) continue;

		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;
		std::string input_name(buffer_filename.str());
		contact_state_reader->read_ContactStatesWOIni(all_frames, input_name);
	}

	cafemol::ClusteringResults clustering_result = analyzer->run<std::array<float, 3>, cafemol::analysis::KMF_INI_PLUSPLUS>(all_frames);
//	cafemol::ClusteringResults clustering_result = analyzer->run<std::array<float, 3>, cafemol::analysis::KMF_INI_DEFAULT>(all_frames);

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Output results to " + output_name);
	// output 

	std::ofstream ofs(output_name, std::ios::out);

	// analysis set-up
	ofs << "<<<< analysis_information" << std::endl;
	ofs << "n_cluster = " << max_cluster_num << std::endl;
	ofs << "max_iteration_num = " << max_iteration_num << std::endl; 
	ofs << "similarity_cutoff = " << cutoff << std::endl;
	ofs << "random_seed = " << n_seed << std::endl; 
	ofs << "bin_width = " << bin_w << std::endl;
	ofs << ">>>>" << std::endl << std::endl;

	// cluster id
	ofs << "<<<< cluster_ids" << std::endl;
	for (std::size_t i_cluster = 0; i_cluster < max_cluster_num; ++i_cluster) {
		const std::vector<std::size_t>& cluster_ids = std::get<0>(clustering_result)[i_cluster];
		ofs << "Cluster: " << i_cluster + 1 << std::endl;
		for (const std::size_t& frame_id : cluster_ids) {
			ofs << frame_id << " ";
		}
		ofs << std::endl;
	}
	ofs << ">>>>" << std::endl << std::endl;

	// intracluster square error
	ofs << "<<<< intra_cluster_se" << std::endl;

	for (std::size_t i_cluster = 0; i_cluster < max_cluster_num; ++i_cluster) {
		const float& intracluster_distance = std::get<1>(clustering_result)[i_cluster];
		ofs << "Cluster" << i_cluster + 1 << ": " << intracluster_distance << std::endl;
	}
	ofs << ">>>>" << std::endl << std::endl;

	ofs << "<<<< similarity_sum" << std::endl;
	float intracluster_ss = 0.0;
	for (std::size_t i_cluster = 0; i_cluster < max_cluster_num; ++i_cluster) {
		intracluster_ss += std::get<1>(clustering_result)[i_cluster];
	}
	ofs << "SS = " << intracluster_ss << std::endl;
	ofs << ">>>>" << std::endl << std::endl;

	// intercluster distance
	ofs << "<<<< inter_cluster_similarity" << std::endl;
	
	ofs << "Cluster ID";
	for (std::size_t i_cluster = 0; i_cluster < max_cluster_num; ++i_cluster) {
		ofs << " " << i_cluster + 1;
	}
	ofs << std::endl;

	const std::vector<std::vector<float>> intercluster_distance = std::get<2>(clustering_result);

	for (std::size_t i_cluster = 0; i_cluster < max_cluster_num; ++i_cluster) {
		ofs << i_cluster + 1 << ":";
		for (std::size_t j_cluster = 0; j_cluster < max_cluster_num; ++j_cluster) {
			ofs << " " << intercluster_distance[i_cluster][j_cluster];
		}
		ofs << std::endl;
	}

	ofs << ">>>>" << std::endl;

	end = std::chrono::system_clock::now();

	double measured_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.);

	std::cout << "Measured Time: " << measured_time << " [s]" << std::endl;

	return 0;
}
