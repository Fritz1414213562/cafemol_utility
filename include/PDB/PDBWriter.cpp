#include"PDBWriter.hpp"


void cafemol::PDBWriter::close_file(std::ifstream& input_file) {
	if (input_file.is_open()) input_file.close();
}


std::string cafemol::PDBWriter::read_LineBeforeCods(const std::string& line) {
	std::string res = read_LineRangeof(line, LINE_BEGIN, BEFORE_COORDINATES);
	return res;
}


std::string cafemol::PDBWriter::read_LineAfterCods(const std::string& line) {
	std::string res;
	if (line.size() > AFTER_COORDINATES) {
		res = read_LineRangeof(line, AFTER_COORDINATES, line.size() - 1);
	}

	else if (line.size() == AFTER_COORDINATES) res = "";

	else if (line.size() < AFTER_COORDINATES) {
		std::cerr << "Error: This row is shorter than 54, so this is not ATOM Column or the column is broken." << std::endl;
		std::exit(1);
	}

	return res;
}


int cafemol::PDBWriter::read_SerialNumof(const std::string& line) {
	std::string buf = read_LineRangeof(line, SERIAL_NUM_BEGIN, SERIAL_NUM_END);
	int res = stoi(buf);
	return res;
}


void cafemol::PDBWriter::replace_CodsInPDBLine(std::ofstream& ofs, const std::string& line, const std::array<float, 3>& vec) {

	std::string line_before_cods = read_LineBeforeCods(line);
	std::string line_after_cods = read_LineAfterCods(line);
	// write string before coordinates columns.
	ofs << line_before_cods; 
	// write coordinates
	for (std::size_t idim = 0; idim < 3; ++idim) {
		// ofs << std::right << std::setw(COORDINATES_COLUMN_SIZE) 
		// << std::setprecision(COORDINATES_COLUMN_SIZE) << vec[idim];
		ofs << std::right << std::setw(8) << std::fixed << std::setprecision(COORDINATES_DECIMAL_SIZE) << vec[idim];
	}
	// write string after coordinates columns.
	ofs << line_after_cods << std::endl;
}


std::vector<std::array<float, 3>> cafemol::PDBWriter::transpose_Mat(const std::array<std::vector<float>, 3>& matrix) {

	// judge the size of matrix
	{
		std::array<int, 3> size_of_each_vector;
		for (std::size_t idim = 0; idim < 3; ++idim) {
			size_of_each_vector[idim] = matrix[idim].size();
		}
		if (size_of_each_vector[0] != size_of_each_vector[1]) {
			std::cerr << "Error: The size of dim0 in the matrix is not consistent with the one of dim1." << std::endl;
			std::exit(1);
		}
		else if (size_of_each_vector[1] != size_of_each_vector[2]) {
			std::cerr << "Error: The size of dim1 in the matrix is not consistent with the one of dim2." << std::endl;
			std::exit(1);
		}
		else {
			std::cout << "The matrix size is the same among 3 dimensions." << std::endl;
		}
	}

	// transpose
	std::vector<std::array<float, 3>> result;

	std::size_t vector_size = matrix[0].size();
	result.resize(vector_size);
	
	for (std::size_t idx = 0; idx < vector_size; ++idx) {
		result[idx] = {matrix[0][idx], matrix[1][idx], matrix[2][idx]};
	}

	return result;
}


void cafemol::PDBWriter::makePDBFileFrom(const std::string& template_name, const std::string& out_name, const std::array<std::vector<float>, 3>& vecs) {

	std::vector<std::array<float, 3>> transposed_vecs = transpose_Mat(vecs);

	std::ifstream ifs(template_name, std::ios::in);
	std::ofstream ofs(out_name, std::ios::out);
	
	// initialize line and index
	std::string line;
	int index = 0;
	int vector_size = transposed_vecs.size();

	while (getline(ifs, line)) {
		if (is_TERRow(line)) {
			ofs << "TER" << std::endl;
			continue;
		}
		else if (!is_AtomRow(line)) continue;

		int atom_ID = read_SerialNumof(line);

		if (atom_ID > vector_size) {
			std::cerr << "The number of atoms in this PDB file is larger than matrix" << std::endl;
			std::exit(1);
		}

		replace_CodsInPDBLine(ofs, line, transposed_vecs[index]);
		++index;
	}

	ofs << "END" << std::endl;
	ifs.close();
	ofs.close();
}
