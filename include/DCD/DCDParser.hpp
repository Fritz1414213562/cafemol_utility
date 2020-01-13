// This class is a parser class for dcd binary file.
//-------------------------------------------------------------------------------
// include guard
#ifndef DCD_PARSER_HPP
#define DCD_PARSER_HPP

//-------------------------------------------------------------------------------
// includes

#include<iostream>
#include<fstream>
#include<vector>
#include<cstdint>
#include<cstring>
#include<array>
#include"../misc/ErrorMessage.hpp"

//-------------------------------------------------------------------------------
// namespace
namespace cafemol {
//-------------------------------------------------------------------------------
// header


class DCDParser {

private :
    
    const std::size_t byte_size = 4;
    std::ifstream inputfile_stream;

protected :

    std::string input_name;
	int frame_num;
	int atom_num;
//	int frame_num;
//	int atom_num;
	cafemol::error_handling::Error_Output error_output = cafemol::error_handling::Error_Output();
    
    template<typename T>
    T read_binary_as(const char *str) {
        T result;
        std::memcpy(std::addressof(result), str, sizeof(T));
        return result;
    }

    void open_file();
    void close_file();
    void load_file(const std::string input_file_name);
	void load_file();
    std::string read_block();
    int read_frame_num(const std::string& first_block);
    int read_atom_num(const std::string& second_block);
    std::vector<float> read_coordinates(const std::string& block, const int atom_num);
    std::array<std::vector<float>, 3> read_xyz(const int atom_num);
	void read_num_frame_and_atom();

public :
    DCDParser();
    DCDParser(const std::string input_file_name);
};

template<typename T, std::size_t N>
inline std::array<T, N> operator-(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    
    std::array<T, N> result;
    for (std::size_t iN = 0; iN < N; ++iN) {
        result[iN] = lhs[iN] - rhs[iN];
    }

    return result;
}


template<typename T, std::size_t N>
inline std::array<T, N> operator*(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    
    std::array<T, N> result;
    for (std::size_t iN = 0; iN < N; ++iN) {
        result[iN] = lhs[iN] * rhs[iN];
    }

    return result;
}


template<typename T, std::size_t N>
inline std::array<T, N>& operator/=(std::array<T, N>& lhs, const T rhs) {
    
    for (std::size_t iN = 0; iN < N; ++iN) {
        lhs[iN] /= rhs;
    }

    return lhs;
}


template<typename T, std::size_t N>
inline std::array<T, N>& operator-=(std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    
    for (std::size_t iN = 0; iN < N; ++iN) {
        lhs[iN] -= rhs[iN];
    }

    return lhs;
}


template<typename T, std::size_t N>
inline std::array<T, N>& operator+=(std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    
    for (std::size_t iN = 0; iN < N; ++iN) {
        lhs[iN] += rhs[iN];
    }

    return lhs;
}


template<typename T, std::size_t N>
inline std::array<T, N> operator/(const std::array<T, N>& lhs, const T rhs) {
    
    std::array<T, N> result;
    for (std::size_t iN = 0; iN < N; ++iN) {
        result[iN] = lhs[iN] / rhs;
    }

    return result;
}


template<typename T, std::size_t N>
inline std::array<T, N> operator*(const std::array<T, N>& lhs, const T rhs) {
    
    std::array<T, N> result;
    for (std::size_t iN = 0; iN < N; ++iN) {
        result[iN] = lhs[iN] * rhs;
    }

    return result;
}
}

#endif /* DCD_PARSER_HPP */
