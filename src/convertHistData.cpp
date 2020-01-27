#include<DensityFunction.hpp>
#include<iostream>
#include<string>
#include<fstream>
#include<array>
#include<vector>
#include<memory>
#include<algorithm>


int main(int argc, char *argv[]) {

	// read command line variants
	std::string input_name = argv[1];
	std::string output_name = argv[2];
	std::string s_bin_num = argv[3];
	bool is_set_data_range;
	float data_range_left_lim;
	float data_range_right_lim;
	if (argc == 6) {
		data_range_left_lim = std::stof(argv[4]);
		data_range_right_lim = std::stof(argv[5]);
		is_set_data_range = true;
	}
	else if (argc == 4) is_set_data_range = false;
	else {
		std::cerr << "Error: too much or less arguments." << std::endl;
		std::exit(1);
	}
	float bin_range = stof(s_bin_num);

	// generate the instance
	std::unique_ptr<cafemol::DensityFunction<float, float>> hist_drawer = std::make_unique<cafemol::DensityFunction<float, float>>();
	// read the input file.
	std::ifstream ifs(input_name, std::ios::in);
	std::vector<float> bp_ids;
	std::string buffer;
	while (std::getline(ifs, buffer)) {
		if (buffer != "nan") {
			float numberof_bp = stof(buffer);
			bp_ids.push_back(numberof_bp);
		}
		else break;
	}
	ifs.close();

	// calculate the density of input vector.
	if (is_set_data_range) hist_drawer->set_DataRange(data_range_left_lim, data_range_right_lim);
	
	std::vector<cafemol::HistogramData<float, float>> hist_xy_data = hist_drawer->make_Density(bp_ids, bin_range);
	const float& vector_size = static_cast<float>(bp_ids.size());

	for (std::size_t idx = 0; idx < hist_xy_data.size(); ++idx) {
		std::get<1>(hist_xy_data[idx]) /= vector_size;
	}

	// make output file
	std::ofstream ofs(output_name, std::ios::out);
//	int hist_x_size = hist_xy_data[0].size();
//	int hist_y_size = hist_xy_data[1].size();
//	if (hist_x_size != hist_y_size) {
//		std::cerr << "Error: The data size of Histgram data in x-axis is not consistent with that in y-axis." << std::endl;
//		std::exit(1);
//	}

	for (std::size_t idx = 0; idx < hist_xy_data.size(); ++idx) {
		ofs << std::get<0>(hist_xy_data[idx]) << " " << std::get<1>(hist_xy_data[idx]) << std::endl;
	}

	ofs.close();

	return 0;
}
