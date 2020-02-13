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
	bool is_set_data_range;
	int data_range_left_lim;
	int data_range_right_lim;
	if (argc == 5) {
		data_range_left_lim = std::stoi(argv[3]);
		data_range_right_lim = std::stoi(argv[4]);
		is_set_data_range = true;
	}
	else if (argc == 3) is_set_data_range = false;
	else {
		std::cerr << "Error: too much or less arguments." << std::endl;
		std::exit(1);
	}

	// generate the instance
	std::unique_ptr<cafemol::DensityFunction<int, std::size_t>> hist_drawer = std::make_unique<cafemol::DensityFunction<int, std::size_t>>();
	// read the input file.
	std::ifstream ifs(input_name, std::ios::in);
	std::vector<int> bp_ids;
	std::string buffer;
	float vector_size = 0;
	while (std::getline(ifs, buffer)) {
		if (buffer != "nan") {
			++vector_size;
			int bp_id = std::stoi(buffer);
			bp_ids.push_back(bp_id);
		}
		else break;
	}
	ifs.close();

	// calculate the density of input vector.
	if (is_set_data_range) hist_drawer->set_DataRange(data_range_left_lim, data_range_right_lim);
	
	std::vector<cafemol::HistogramData<int, std::size_t>> hist_xy_data = hist_drawer->make_Density(bp_ids);

	// make output file
	std::ofstream ofs(output_name, std::ios::out);

	for (std::size_t idx = 0; idx < hist_xy_data.size(); ++idx) {
		ofs << std::get<0>(hist_xy_data[idx]) << " " 
			<< static_cast<float>(std::get<1>(hist_xy_data[idx])) / vector_size << std::endl;
	}

	ofs.close();

	return 0;
}
