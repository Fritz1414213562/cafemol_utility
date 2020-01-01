#ifndef DCD_WRITER_HPP
#define DCD_WRITER_HPP
#include"DCDParser.hpp"
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<tuple>


namespace cafemol {

class DCDWriter : DCDParser {

public:
	
	DCDWriter(const std::string& output_file_name);
	DCDWriter(const std::string& output_file_name, const std::string& input_file_name);


	// test
	void copy_DCD(int& istart, int& nstep_save, int& number_of_steps, int& number_of_units, float& delta, int& version);

private:

	enum BLOCK_SIZES {
		// block size
		first_block_size = 84,
		second_block_size = 324,
		third_block_size = 4,

		// block unit
	};

	char block_char = '=';
	char block_space = ' ';
	char CORD_C = 'C';
	char CORD_O = 'O';
	char CORD_R = 'R';
	char CORD_D = 'D';
	int i_zero = 0;
	int block_unit = 4;

	std::string output_name;


	template<typename T>
	void write_Binary(T& value, std::ofstream& ofs) {
		ofs.write(reinterpret_cast<char*>(&value), sizeof(T));
	}


	template<std::size_t iterative_num = 0, typename Tuple_Type>
	void iterate_Binary(Tuple_Type& Tuple, std::ofstream& ofs) {
		if constexpr ((iterative_num >= 0) && (iterative_num < std::tuple_size<Tuple_Type>::value)) {
			ofs.write(reinterpret_cast<char*>(&std::get<iterative_num>(Tuple)), sizeof(std::get<iterative_num>(Tuple)));
			iterate_Binary<iterative_num + 1>(Tuple, ofs);
		}
	}


//	template<BLOCK_SIZES block_number, typename ...Types>
	template<typename ...Types>
	void write_Block(std::ofstream& ofs, Types& ...arguments) {

		std::tuple<Types...> args_tuple = std::make_tuple(arguments...);
		// calculate the total type size of tuple elements
//		constexpr int tuple_elements_size = count_TupleElementsSize(args_tuple);
		int tuple_elements_size = count_TupleElementsSize(args_tuple);
		// check the block_size
//		static_assert(block_number == tuple_elements_size, "Block size is not consistent with tuple elements size");

		std::cout << tuple_elements_size << std::endl;
		// block_size
		write_Binary(tuple_elements_size, ofs);
		// write body of the block
		iterate_Binary(args_tuple, ofs);
		// check-block size
		write_Binary(tuple_elements_size, ofs);
	}


	// compile-time methods

	template<std::size_t iterative_num = 0, typename Tuple_Type>
	constexpr int count_TupleElementsSize(Tuple_Type& Tuple) {

		static_assert((0 <= iterative_num) || (iterative_num < std::tuple_size<Tuple_Type>::value), "Invalid Index, out of tuple range");

		if constexpr (iterative_num == std::tuple_size<Tuple_Type>::value - 1) {
			return sizeof(std::get<iterative_num>(Tuple));
		}
		else if constexpr (iterative_num < std::tuple_size<Tuple_Type>::value - 1) {
			return sizeof(std::get<iterative_num>(Tuple)) + count_TupleElementsSize<iterative_num + 1>(Tuple);
		}
	}


	// for writing each block
	void write_1stBlock(std::ofstream& ofs, int& frame_number, int& istart, int& nstep_save, int& number_of_steps, int& number_of_units, float& delta, int& version) {
//	void write_1stBlock(std::ofstream& ofs, const int& frame_number, const int& istart, const int& nstep_save, const int& number_of_steps, const int& number_of_units, const float& delta, const int& version);
		write_Block(ofs, CORD_C, CORD_O, CORD_R, CORD_D, frame_number, istart, nstep_save, number_of_steps, number_of_units, i_zero, i_zero, i_zero, i_zero, delta, i_zero, i_zero, i_zero, i_zero, i_zero, i_zero, i_zero, i_zero, i_zero, version);
	}

	void write_2ndBlock(std::ofstream& ofs) {
//	void write_2ndBlock(std::ofstream& ofs);
		write_Block(ofs, block_unit, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char, block_char);
	}

	void write_3rdBlock(std::ofstream& ofs, int& atom_number) {
//	void write_3rdBlock(std::ofstream& ofs,  int& atom_number);
		write_Block(ofs, atom_number);
	}

//
	void write_Body(std::ofstream& ofs, std::array<std::vector<float>, 3>& xyz);

};
}

#endif /* DCD_WRITER_HPP */
