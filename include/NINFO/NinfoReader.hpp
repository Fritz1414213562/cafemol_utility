#ifndef NINFO_READER_HPP
#define NINFO_READER_HPP

#include"NinfoParser.hpp"
#include<array>
#include<vector>
#include<iostream>

namespace cafemol {

class NinfoReader : public NinfoParser {

public:
	NinfoReader() = default;
	NinfoReader(const std::string& input_file_name);
	~NinfoReader() = default;

	enum_types::ninfo_column convert_LineKind2Id(const std::string& line_kind);

	ninfo_data_type::container_bond read_NinfoBond();
	ninfo_data_type::container_angl read_NinfoAngl();
	ninfo_data_type::container_aicg13 read_NinfoAicg13();
	ninfo_data_type::container_dihd read_NinfoDihd();
	ninfo_data_type::container_aicgdih read_NinfoAichdih();
	ninfo_data_type::container_contact read_NinfoContact();
	ninfo_data_type::container_pdns read_NinfoPdns();


private:
	void judge_FileOpen();

	template<enum_types::ninfo_column line_kind>
	std::vector<typename ninfo_data_type::nativeinfo_tuple<line_kind>> read_NinfoData() {
//	std::vector<ninfo_data_type::nativeinfo_tuple<line_kind>> read_NinfoData() {
		judge_FileOpen();

		std::string line;
		std::vector<typename ninfo_data_type::nativeinfo_tuple<line_kind>> data;
		while (std::getline(input_file, line)) {
			std::string head_line = get_Headof(line);
			if (head_line.empty()) continue;
			else if (convert_LineKind2Id(head_line) != line_kind) continue;
			else if (convert_LineKind2Id(head_line) == line_kind) {
				ninfo_data_type::nativeinfo_tuple<line_kind> res = read_Line<line_kind>(line);
				data.push_back(res);
			}
		}
		// reload file
		load_file();
		
		return data;
	}

};

}

#endif /* NINFO_READER_HPP */
