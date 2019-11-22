#include"NinfoReader.hpp"


cafemol::NinfoReader::NinfoReader(const std::string& input_file_name) : NinfoParser(input_file_name) {}


cafemol::enum_types::ninfo_column cafemol::NinfoReader::convert_LineKind2Id(const std::string& line_kind) {

	if (line_kind == "bond") return enum_types::bond;
	else if (line_kind == "angl") return enum_types::angl;
	else if (line_kind == "aicg13") return enum_types::aicg13;
	else if (line_kind == "dihd") return enum_types::dihd;
	else if (line_kind == "aicgdih") return enum_types::aicgdih;
	else if (line_kind == "contact") return enum_types::contact;
	else if (line_kind == "pdns") return enum_types::pdns;
	else return enum_types::other;

}


void cafemol::NinfoReader::judge_FileOpen() {
	if (!input_file.is_open()) {
		std::cerr << "Error: File is not opened" << std::endl;
		std::exit(1);
	}
}


//std::vector<cafemol::ninfo_data_type::nativeinfo_tuple<cafemol::enum_types::bond>>
//std::vector<cafemol::ninfo_data_type::bond_tuple>
cafemol::ninfo_data_type::container_bond cafemol::NinfoReader::read_NinfoBond() {
	return read_NinfoData<enum_types::bond>();
}


cafemol::ninfo_data_type::container_angl cafemol::NinfoReader::read_NinfoAngl() {
	return read_NinfoData<enum_types::angl>();
}


cafemol::ninfo_data_type::container_aicg13 cafemol::NinfoReader::read_NinfoAicg13() {
	return read_NinfoData<enum_types::aicg13>();
}


cafemol::ninfo_data_type::container_dihd cafemol::NinfoReader::read_NinfoDihd() {
	return read_NinfoData<enum_types::dihd>();
}


cafemol::ninfo_data_type::container_aicgdih cafemol::NinfoReader::read_NinfoAichdih() {
	return read_NinfoData<enum_types::aicgdih>();
}


cafemol::ninfo_data_type::container_contact cafemol::NinfoReader::read_NinfoContact() {
	return read_NinfoData<enum_types::contact>();
}


cafemol::ninfo_data_type::container_pdns cafemol::NinfoReader::read_NinfoPdns() {
	return read_NinfoData<enum_types::pdns>();
}

