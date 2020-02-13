#ifndef MATRIX_FILE_READER_HPP
#define MATRIX_FILE_READER_HPP
#include<OtherFormat/MatrixFileParser.hpp>
#include<string>
#include<array>
#include<Eigen/Core>

namespace cafemol {

class MatrixFileReader : public MatrixFileParser {

public:

	MatrixFileReader() : MatrixFileParser() {}
	MatrixFileReader(const std::string& inputfile_name) : MatrixFileParser(inputfile_name) {}
	~MatrixFileReader() = default;


	Eigen::MatrixXd read_Matrix() {
		open_File();

		std::array<int, 2> mat_shape = read_MatShape();
		Eigen::MatrixXd result(mat_shape[0], mat_shape[1]);
		std::string s_mat_data = read_Block();
		close_File();

		for (int i_datum = 0; i_datum < result.size(); ++i_datum) {
			std::size_t pos_of_datum = i_datum * sizeof(double);
			result(i_datum) = read_Binary_as<double>(&s_mat_data.at(pos_of_datum));
		}


		return result;
	}



	Eigen::MatrixXd read_Matrix(const std::string& input_name) {
		open_File(input_name);

		std::array<int, 2> mat_shape = read_MatShape();
		const int& mat_size = mat_shape[0] * mat_shape[1];
		std::cout << "Row Size -> " << mat_shape[0] << std::endl;
		std::cout << "Column Size -> " << mat_shape[1] << std::endl << std::endl;

		Eigen::MatrixXd result(mat_shape[0], mat_shape[1]);
		std::string s_mat_data = read_Block();
		close_File();
		std::cout << input_name << " closed" << std::endl;

		for (int i_datum = 0; i_datum < mat_size; ++i_datum) {
			std::size_t pos_of_datum = i_datum * sizeof(double);
			result(i_datum) = read_Binary_as<double>(&s_mat_data.at(pos_of_datum));
		}


		return result;
	}

};

}

#endif /* MATRIX_FILE_READER_HPP */
