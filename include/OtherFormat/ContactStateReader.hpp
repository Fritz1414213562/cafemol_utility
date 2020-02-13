#ifndef CONTACT_STATE_READER_HPP
#define CONTACT_STATE_READER_HPP
#include<OtherFormat/ContactStateParser.hpp>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

namespace cafemol {

class ContactStateReader : public ContactStateParser {

public:

	ContactStateReader() : ContactStateParser() {}
	
	ContactStateReader(const std::string& input_filename) : ContactStateParser(input_filename) {}
	~ContactStateReader() = default;


	void read_ContactStatesWOIni(std::vector<std::vector<std::array<int, 2>>>& all_data, const std::string& file_name) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<std::array<int,2>> data_on_frame;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && !data_on_frame.empty()) {
				all_data.push_back(data_on_frame);
				data_on_frame.clear();
				continue;
			}

			if (is_reading_1stFrame) continue;

			data_on_frame.push_back({std::stoi(buffer_words[0]),
									 std::stoi(buffer_words[1])});
		}

		if (!data_on_frame.empty()) {
			all_data.push_back(data_on_frame);
		}

		close_File();
	}


	void read_ContactStatesWOIni(std::vector<std::vector<float>>& all_data, const std::array<std::size_t, 2>& data_frame_size, const std::string& file_name) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<float> data_on_frame(data_frame_size[0] * data_frame_size[1], 0.0);

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}

			if (is_reading_1stFrame) continue;
		//	else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
			if ((buffer_words[0] == "Frame") && !is_reading_1stFrame) {
				all_data.push_back(data_on_frame);
				data_on_frame.clear();
				data_on_frame.resize(data_frame_size[0] * data_frame_size[1], 0.0);
				continue;
			}


			const std::size_t& data_coord_x = std::stoi(buffer_words[0]) - 1;
			const std::size_t& data_coord_y = std::stoi(buffer_words[1]) - 1;

			data_on_frame[data_coord_y * data_frame_size[0] + data_coord_x] = 1.0;

		}

		if (!data_on_frame.empty()) {
			all_data.push_back(data_on_frame);
		}

		close_File();
	}






	void read_ContactStatesWOIni(std::vector<std::vector<float>>& all_data, const std::array<std::size_t, 2>& data_frame_size, const std::string& file_name, const std::size_t& interval) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<float> data_on_frame(data_frame_size[0] * data_frame_size[1], 0.0);
	//	std::size_t frame_num = 0;
		bool is_on_frame = false;
		std::size_t read_num = 0;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}

			if ((buffer_words[0] == "Frame") && (std::stoi(buffer_words[1]) % interval == 0)) {
				is_on_frame = true;
				continue;
			}
//			else if ((buffer_words[0] == "Frame") && (std::stoi(buffer_words[1]) % interval != 0)) {
//				is_on_frame = false;
//				continue;
//			}


			if (is_reading_1stFrame) continue;
			else if (!is_on_frame) continue;
		//	else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
			if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && is_on_frame) {
				all_data.push_back(data_on_frame);
				data_on_frame.clear();
				data_on_frame.resize(data_frame_size[0] * data_frame_size[1], 0.0);
				++read_num;
				is_on_frame = false;
				continue;
			}


			const std::size_t& data_coord_x = std::stoi(buffer_words[0]) - 1;
			const std::size_t& data_coord_y = std::stoi(buffer_words[1]) - 1;

			data_on_frame[data_coord_y * data_frame_size[0] + data_coord_x] = 1.0;

		}

		if (!data_on_frame.empty()) {
			all_data.push_back(data_on_frame);
			++read_num;
		}

	//	std::cout << read_num << std::endl;

		close_File();
	}



	void read_ContactStatesWOIni(std::vector<std::vector<std::array<int, 2>>>& all_data, const std::string& file_name, const std::size_t& interval) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<std::array<int, 2>> data_on_frame;
	//	std::size_t frame_num = 0;
		bool is_on_frame = false;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}

			if ((buffer_words[0] == "Frame") && (std::stoi(buffer_words[1]) % interval == 0)) {
				is_on_frame = true;
				continue;
			}
//			else if ((buffer_words[0] == "Frame") && (std::stoi(buffer_words[1]) % interval != 0)) {
//				is_on_frame = false;
//				continue;
//			}


			if (is_reading_1stFrame) continue;
			else if (!is_on_frame) continue;
		//	else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
			if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && is_on_frame) {
				all_data.push_back(data_on_frame);
				data_on_frame.clear();
				is_on_frame = false;
				continue;
			}


			const int& data_coord_x = std::stoi(buffer_words[0]);
			const int& data_coord_y = std::stoi(buffer_words[1]);

			data_on_frame.push_back({data_coord_x, data_coord_y});

		}

		if (!data_on_frame.empty()) {
			all_data.push_back(data_on_frame);
		}

	//	std::cout << read_num << std::endl;

		close_File();
	}









	std::vector<std::vector<std::array<int, 2>>> read_ContactStatesWOIni(const std::string& file_name) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<std::vector<std::array<int, 2>>> result;
		std::vector<std::array<int,2>> data_on_frame;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && !data_on_frame.empty()) {
				result.push_back(data_on_frame);
				data_on_frame.clear();
				continue;
			}

			if (is_reading_1stFrame) continue;

			data_on_frame.push_back({std::stoi(buffer_words[0]),
									 std::stoi(buffer_words[1])});
		}

		if (!data_on_frame.empty()) {
			result.push_back(data_on_frame);
		}

		close_File();

		return result;
	}



	template<typename Size3_Vector>
	void read_ContactStatesWOIni(std::vector<std::vector<Size3_Vector>>& all_data, const std::string& file_name) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<Size3_Vector> data_on_frame;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && data_on_frame.empty()) continue;
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && !data_on_frame.empty()) {
				all_data.push_back(data_on_frame);
				data_on_frame.clear();
				continue;
			}

			if (is_reading_1stFrame) continue;

			data_on_frame.push_back({std::stof(buffer_words[0]),
									 std::stof(buffer_words[1]),
									 1.0});
		}

		if (!data_on_frame.empty()) {
			all_data.push_back(data_on_frame);
		}

		close_File();
	}


	void read_ContactStatesFrameSizeWOIni(std::vector<std::size_t>& all_data_frame_size, const std::string& file_name) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::size_t data_frame_size = 0;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check", file_name);
			if ((buffer_words[0] == "Frame") && (buffer_words[1] == "0")) {
				is_reading_1stFrame = true;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && is_reading_1stFrame) {
				is_reading_1stFrame = false;
				continue;
			}
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && (data_frame_size == 0)) continue;
			else if ((buffer_words[0] == "Frame") && !is_reading_1stFrame && (data_frame_size > 0)) {
				all_data_frame_size.push_back(data_frame_size);
				data_frame_size = 0;
				continue;
			}

			if (is_reading_1stFrame) continue;

			++data_frame_size;
		}

		if (data_frame_size > 0) {
			all_data_frame_size.push_back(data_frame_size);
		}

		close_File();
	}




};

}


#endif /* CONTACT_STATE_READER_HPP */
