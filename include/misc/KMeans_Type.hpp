#ifndef KMEANS_TYPE_HPP
#define KMEANS_TYPE_HPP


namespace cafemol::analysis {


template<typename realT>
struct FrameData {

public:

	FrameData(const std::size_t& row, const std::size_t& col) : row_size(row), column_size(col) {
		data.resize(row * col);
	}
	~FrameData() = default;

	const FrameData& operator=(const std::vector<std::vector<realT>>& input_data) noexcept {
		data = input_data;
	}

	FrameData operator=(std::vector<std::vector<realT>> input_data) noexcept {
		data = input_data;
	}

	const std::vector<realT>& operator[](const std::size_t& input_index) {
		return data[input_index];
	}

	std::vector<realT> operator[](const std::size_t& input_index) {
		return data[input_index];
	}

	const realT& operator[](const std::size_t& i_row, const std::size_t& i_col) noexcept {
		return data[i_row][i_col];
	}

	realT operator[](const std::size_t& i_row, const std::size_t& i_col) noexcept {
		return data[i_row][i_col];
	}

	const std::size_t& rows() {return row_size;}
	std::size_t rows() {return row_size;}
	const std::size_t& cols() {return column_size;}
	std::size_t cols() {return column_size;}

private:
	std::vector<std::vector<realT>> data;

	const std::size_t row_size;
	const std::size_t column_size;


};

}


#endif /* KMEANS_TYPE_HPP */
