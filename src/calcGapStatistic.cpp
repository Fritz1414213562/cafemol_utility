#include<OtherFormat/ClusteringResultReader.hpp>
#include<ErrorMessage.hpp>
#include<StandardOutput.hpp>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<memory>
#include<array>
#include<vector>
#include<string>


int main(int argc, char *argv[]) {

	const int& max_argc = 8;
	cafemol::error_handling::Error_Output eout = cafemol::error_handling::Error_Output();
	if (argc != max_argc) eout("too much or less arguments");
	cafemol::output_handling::Standard_Output sout = cafemol::output_handling::Standard_Output();

	const std::string& input_prefix = argv[1];
	const std::string& input_suffix = argv[2];
	const std::string& dummy_prefix = argv[3];
	const std::string& dummy_suffix = argv[4];
	const std::string& s_file_num1 = argv[5];
	const std::string& s_file_num2 = argv[6];
	const std::size_t& i_file_num1 = std::stoi(argv[5]);
	const std::size_t& i_file_num2 = std::stoi(argv[6]);
	const std::array<std::size_t, 2>& filenum_range = {std::min(i_file_num1,
																i_file_num2),
													   std::max(i_file_num1,
																i_file_num2)};
	const std::string& output_name = argv[7];

	const std::size_t& BLOCK_SIZE = 90;
	const std::size_t& file_num_digit = std::max(s_file_num1.size(), s_file_num2.size());

	sout.output_HyphenBlock("", BLOCK_SIZE);
	sout("Reading Clustering Result Data");

	std::ofstream ofs(output_name, std::ios::out);

	for (std::size_t file_index = filenum_range[0]; file_index <= filenum_range[1]; ++file_index) {
		std::stringstream buffer_filename;
		buffer_filename << input_prefix
						<< std::setw(file_num_digit)
						<< std::setfill('0')
						<< file_index
						<< input_suffix;
		std::string input_name(buffer_filename.str());

		std::stringstream buffer_dummyfilename;
		buffer_dummyfilename << dummy_prefix
							 << std::setw(file_num_digit)
							 << std::setfill('0')
							 << file_index
							 << dummy_suffix;

		std::string dummy_name(buffer_dummyfilename.str());

		std::unique_ptr<cafemol::ClusteringResultReader> inputfile_reader = std::make_unique<cafemol::ClusteringResultReader>(input_name);
		std::unique_ptr<cafemol::ClusteringResultReader> dummyfile_reader = std::make_unique<cafemol::ClusteringResultReader>(dummy_name);

		const float& input_similarity_sum = inputfile_reader->read_SimilaritySum();
		const float& dummy_similarity_sum = dummyfile_reader->read_SimilaritySum();

	//	const float gap_statistic = std::log10(dummy_similarity_sum / input_similarity_sum);
		const float gap_statistic = std::log10(dummy_similarity_sum) -
									std::log10(input_similarity_sum);
	//	const float gap_statistic = log10(dummy_similarity_sum / input_similarity_sum);
		ofs << file_index << " " << gap_statistic << std::endl;

	}

	return 0;
}
