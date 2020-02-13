#include<OtherFormat/ContactStateReader.hpp>
#include<analysis/KMeans>
#include<IO/ErrorMessage.hpp>
#include<IO/StandardOutput.hpp>
#include<string>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<memory>
#include<array>
#include<vector>
#include<chrono>


int main(int argc, char *argv[]) {

//	const int& max_argc = 11;
	const int& max_argc = 10;
	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	CafeInLess::IO::Standard_Output sout = CafeInLess::IO::Standard_Output();

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
	const std::size_t data_frame_x = std::stoi(argv[7]);
	const std::size_t data_frame_y = std::stoi(argv[8]);
	const std::array<std::size_t, 2> data_frame_shape = {data_frame_x, data_frame_y};


	const std::size_t n_seed = std::stoi(argv[9]);
//	const std::size_t interval_length = std::stoi(argv[10]);

	const std::size_t BLOCK_SIZE = 90;

	// check a digit of file number
	const std::string s_file_num = argv[2];
	const std::size_t& file_num_digit = s_file_num.size();


	// check Skip Number
	const std::string& skip_filename = "SkipNum";
	std::ifstream skip_file(skip_filename, std::ios::in);
	std::vector<std::size_t> skip_number;
	if (skip_file.is_open()) {
		sout("The file, SkipNum has read.");

		std::string buffer;
		while (std::getline(skip_file, buffer)) {
			skip_number.push_back(std::stoi(buffer));
		}
	}
	skip_file.close();


// body

	std::unique_ptr<cafemol::ContactStateReader> contact_state_reader = std::make_unique<cafemol::ContactStateReader>();
	std::unique_ptr<CafeInLess::analysis::KMeans<CafeInLess::analysis::KMEANS_DIST_L2>> k_means_device = std::make_unique<CafeInLess::analysis::KMeans<CafeInLess::analysis::KMEANS_DIST_L2>>(max_cluster_num, max_iteration_num, n_seed);

//	FrameData<float> frames;
//	{
	std::vector<std::vector<float>> all_frames;

	sout[BLOCK_SIZE];

	sout("Calculation starts.", "");

	sout("Reading Contact State Data");


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
	//	contact_state_reader->read_ContactStatesWOIni(all_frames, data_frame_shape, input_name, interval_length);
		contact_state_reader->read_ContactStatesWOIni(all_frames, data_frame_shape, input_name);
	}
//		frames = all_frames
//	}

//	const std::vector<std::vector<std::size_t>>& clustering_result = k_means_device->run(frames);
	const std::vector<std::vector<std::size_t>>& clustering_result = k_means_device->run(all_frames);
	sout("", "Calculation ends.");

	sout[BLOCK_SIZE];
	sout("Output results to " + output_name);
	// output 

	std::ofstream ofs(output_name, std::ios::out);

	// analysis set-up
	ofs << "<<<< analysis_information" << std::endl;
	ofs << "n_cluster = " << clustering_result.size() << std::endl;
	std::size_t c_id = 1;
	for (const std::vector<std::size_t>& clustering_ids : clustering_result) {
		ofs << "cluster_id" << c_id << " " << clustering_ids.size() << std::endl;
		++c_id;
	}

	ofs << "max_iteration_num = " << max_iteration_num << std::endl; 
	ofs << "random_seed = " << n_seed << std::endl; 
	ofs << ">>>>" << std::endl << std::endl;

	// cluster id
	ofs << "<<<< cluster_ids" << std::endl;
	for (std::size_t i_cluster = 0; i_cluster < clustering_result.size(); ++i_cluster) {
		const std::vector<std::size_t>& cluster_ids = clustering_result[i_cluster];
		ofs << "Cluster: " << i_cluster + 1 << std::endl;
		for (const std::size_t& frame_id : cluster_ids) {
			ofs << frame_id << " ";
		}
		ofs << std::endl;
	}
	ofs << ">>>>" << std::endl << std::endl;


	end = std::chrono::system_clock::now();

	double measured_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.);

	std::cout << "Measured Time: " << measured_time << " [s]" << std::endl;



	return 0;
}
