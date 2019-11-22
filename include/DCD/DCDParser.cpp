#include"DCDParser.hpp"

//-------------------------------------------------------------------------------
// body

cafemol::DCDParser::DCDParser() {
    close_file();    
}


cafemol::DCDParser::DCDParser(const std::string input_file_name) {
    input_name = input_file_name;

    if (inputfile_stream.is_open()) {
        inputfile_stream.close();
    }
}


void cafemol::DCDParser::open_file() {

    inputfile_stream.open(input_name, std::ios::in | std::ios::binary);

    if (!inputfile_stream) {
        std::cerr << "FileOpenError: File does not exist." << std::endl;
        std::exit(1);
    }
}


void cafemol::DCDParser::close_file() {
    if (inputfile_stream.is_open()) {
        inputfile_stream.close();
     }
}


void cafemol::DCDParser::load_file(const std::string input_file_name) {

	close_file();
    input_name = input_file_name;
}


void cafemol::DCDParser::load_file() {
	close_file();
}


std::string cafemol::DCDParser::read_block() {
    // read first block size in the line and the data 
    std::int32_t block_size;
    std::vector<char> buffer;
    //char* buffer;
    //std::string buffer;
    //std::exit(1);
    inputfile_stream.read(reinterpret_cast<char*>(&block_size), 4);
    buffer.resize(block_size);
    inputfile_stream.read(buffer.data(), block_size);

    // read last block size in the line
    std::int32_t check_block_size;
    inputfile_stream.read(reinterpret_cast<char*>(&check_block_size), 4);

    if (block_size != check_block_size) {
        std::cerr << "Error: Left Blocksize is not consistent to Right one." << std::endl;
		std::cerr << "Left -> " << block_size << std::endl;
		std::cerr << "Right -> " << check_block_size << std::endl;
        std::exit(1);
    }

    std::string result(buffer.begin(), buffer.end());
    //std::string result = buffer;

    return result; 
}


int cafemol::DCDParser::read_frame_num(const std::string& first_block) {
    int frame_num = read_binary_as<int>(&first_block.at(byte_size));
    return frame_num;
}


int cafemol::DCDParser::read_atom_num(const std::string& second_block) {
    int atom_num = read_binary_as<int>(&second_block.at(0));
    return atom_num;
}


std::vector<float> cafemol::DCDParser::read_coordinates(const std::string& block, const int atom_num) {
    std::vector<float> coordinates(atom_num);
    for (std::size_t atom_i = 0; atom_i < atom_num; ++atom_i) {
        std::size_t pos_of_atom = atom_i * sizeof(float);
        coordinates[atom_i] = read_binary_as<float>(&block.at(pos_of_atom));
    }
    return coordinates;
}


std::array<std::vector<float>, 3> cafemol::DCDParser::read_xyz(const int atom_num) {
    std::string x_block = read_block();
    std::string y_block = read_block();
    std::string z_block = read_block();
    std::array<std::vector<float>, 3> xyzs;
    xyzs[0] = read_coordinates(x_block, atom_num);
    xyzs[1] = read_coordinates(y_block, atom_num);
    xyzs[2] = read_coordinates(z_block, atom_num);

    return xyzs;
}


void cafemol::DCDParser::read_num_frame_and_atom() {

    open_file();
    // read 1st block in header and get the number of frames
    std::string first_block = read_block();
    frame_num = read_frame_num(first_block);
    // read 2nd block in header
    read_block();
    // read 3rd block in header and get the number of atoms
    std::string third_block = read_block();
    atom_num = read_atom_num(third_block);
    // option: output to standard
    std::cout << "frame_num " << frame_num << std::endl;
    std::cout << "atom_num " << atom_num << std::endl;
}
