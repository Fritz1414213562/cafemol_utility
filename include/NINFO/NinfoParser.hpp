#ifndef NINFO_PARSER_HPP
#define NINFO_PARSER_HPP
#include<iostream>
#include<string>
#include<fstream>
#include<array>
#include<vector>
#include<tuple>

namespace cafemol {

namespace enum_types {

enum ninfo_column {
	bond, angl, aicg13, dihd, aicgdih, contact, pdns, other,
};

}

namespace ninfo_data_type {

	template<unsigned int unit_pair, unsigned int mp_pairs, unsigned int paramer_sets>
	using ninfo_line_data = std::tuple<int, std::array<int, unit_pair>, std::array<int, mp_pairs>, std::array<float, paramer_sets>>;


	template<enum_types::ninfo_column line_kind>
	struct nativeinfo_tuple;

	template<>
	struct nativeinfo_tuple<enum_types::bond> {
		std::string line_kind_name;
		ninfo_line_data<2, 4, 4> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::angl> {
		std::string line_kind_name;
		ninfo_line_data<2, 6, 4> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::aicg13> {
		std::string line_kind_name;
		ninfo_line_data<2, 6, 5> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::dihd> {
		std::string line_kind_name;
		ninfo_line_data<2, 8, 5> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::aicgdih> {
		std::string line_kind_name;
		ninfo_line_data<2, 8, 5> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::contact> {
		std::string line_kind_name;
		ninfo_line_data<2, 4, 4> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::pdns> {
		std::string line_kind_name;
		ninfo_line_data<1, 2, 4> line_data;
	};

	template<>
	struct nativeinfo_tuple<enum_types::other> {
		std::string line;
	};


	using bond_tuple = nativeinfo_tuple<enum_types::bond>;
	using angl_tuple = nativeinfo_tuple<enum_types::angl>;
	using aicg13_tuple = nativeinfo_tuple<enum_types::aicg13>;
	using dihd_tuple = nativeinfo_tuple<enum_types::dihd>;
	using aicgdih_tuple = nativeinfo_tuple<enum_types::aicgdih>;
	using contact_tuple = nativeinfo_tuple<enum_types::contact>;
	using pdns_tuple = nativeinfo_tuple<enum_types::pdns>;


//	using container_bond = std::vector<nativeinfo_tuple<enum_types::bond>::line_data>;
	using container_bond = std::vector<bond_tuple>;
//	using container_angl = std::vector<nativeinfo_tuple<enum_types::angl>::line_data>;
	using container_angl = std::vector<angl_tuple>;
//	using container_aicg13 = std::vector<nativeinfo_tuple<enum_types::aicg13>::line_data>;
	using container_aicg13 = std::vector<aicg13_tuple>;
//	using container_dihd = std::vector<nativeinfo_tuple<enum_types::dihd>::line_data>;
	using container_dihd = std::vector<dihd_tuple>;
//	using container_aicgdih = std::vector<nativeinfo_tuple<enum_types::aicgdih>::line_data>;
	using container_aicgdih = std::vector<aicgdih_tuple>;
//	using container_contact = std::vector<nativeinfo_tuple<enum_types::contact>::line_data>;
	using container_contact = std::vector<contact_tuple>;
//	using container_pdns = std::vector<nativeinfo_tuple<enum_types::pdns>::line_data>;
	using container_pdns = std::vector<pdns_tuple>;
}


class NinfoParser {

public:
	NinfoParser() = default;
	NinfoParser(const std::string& input_file_name);
	~NinfoParser() = default;

	void load_file(const std::string& file_name);

protected:
	
	std::ifstream input_file;
	std::string input_name;

	// open / close file
	void load_file();

	void close_file();

	std::string get_Headof(const std::string& line);

	template<enum_types::ninfo_column line_kind>
	ninfo_data_type::nativeinfo_tuple<line_kind> read_Line(const std::string& line);


private:

	bool is_Block(const std::string& line);

	std::vector<std::string> split_Line(const std::string& line, const char& delimiter);

	const char ninfo_delimiter = ' ';

};
}

#endif /* NINFO_PARSER_HPP */
