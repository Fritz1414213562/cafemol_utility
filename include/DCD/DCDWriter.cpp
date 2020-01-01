#include"DCDWriter.hpp"


cafemol::DCDWriter::DCDWriter(const std::string& output_file_name) : DCDParser() {
	output_name = output_file_name;
}


cafemol::DCDWriter::DCDWriter(const std::string& output_file_name, const std::string& input_file_name) : DCDParser(input_file_name) {
	output_name = output_file_name;
}


// test
void cafemol::DCDWriter::copy_DCD(int& istart, int& nstep_save, int& number_of_steps, int& number_of_units, float& delta, int& version) {

	if (output_name.empty()) {
		std::cerr << "Error: Output file has not been loaded" << std::endl;
		std::exit(1);
	}
	std::ofstream ofs(output_name, std::ios::out | std::ios::binary);

	// read the headers, frame number, and atom number of input file.
	read_num_frame_and_atom();
	// write the headers
	write_1stBlock(ofs, frame_num, istart, nstep_save, number_of_steps, number_of_units, delta, version);
	write_2ndBlock(ofs);
	write_3rdBlock(ofs, atom_num); 
	
	for (std::size_t iframe = 0; iframe < frame_num; ++iframe) {
		std::array<std::vector<float>, 3> xyz = read_xyz(atom_num);
		write_Body(ofs, xyz);
	}
	load_file();
	ofs.close();
}



//void cafemol::DCDWriter::write_1stBlock(std::ofstream& ofs, const int& frame_number, const int& istart, const int& nstep_save, const int& number_of_steps, const int& number_of_units, const float& delta, const int& version) {
////	write_Block<first_block_size>(ofs, 'C', 'O', 'R', 'D', frame_number, istart, nstep_save, number_of_steps, number_of_units, 0, 0, 0, 0, delta, 0, 0, 0, 0, 0, 0, 0, 0, 0, version);
//	write_Block(ofs, 'C', 'O', 'R', 'D', frame_number, istart, nstep_save, number_of_steps, number_of_units, 0, 0, 0, 0, delta, 0, 0, 0, 0, 0, 0, 0, 0, 0, version);
//}


//void cafemol::DCDWriter::write_2ndBlock(std::ofstream& ofs) {
//	// crazy
////	write_Block<second_block_size>(ofs, 4, '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=');
//	write_Block(ofs, 4, '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=', '=');
//}


//void cafemol::DCDWriter::write_3rdBlock(std::ofstream& ofs, const int& atom_number) {
//	write_Block(ofs, atom_number);
////	write_Block<third_block_size>(ofs, atom_number);
//}


void cafemol::DCDWriter::write_Body(std::ofstream& ofs,  std::array<std::vector<float>, 3>& xyz) {
	for (std::size_t idim = 0; idim < 3; ++idim) {
		int atom_number = xyz[idim].size();
		int block_size = atom_number * sizeof(int);
		// block size
		write_Binary(block_size, ofs);
		// write body of the block
		for (int idx = 0; idx < atom_number; ++idx) {
			write_Binary(xyz[idim][idx], ofs);
		}
		// check-block size
		write_Binary(block_size, ofs);
	}
}
