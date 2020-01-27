#ifndef CONTACT_STATE_READER_HPP
#define CONTACT_STATE_READER_HPP
#include<OtherFormat/ContactStateParser.hpp>
#include<opencv2/core.hpp>
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

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check " + file_name);
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



	template<typename Size3_Vector>
	void read_ContactStatesWOIni(std::vector<std::vector<Size3_Vector>>& all_data, const std::string& file_name) {
		open_File(file_name);
		std::string buffer;
		bool is_reading_1stFrame = false;
		std::vector<Size3_Vector> data_on_frame;

		while(std::getline(input_file, buffer)) {
			const std::vector<std::string> buffer_words = split_ConLine(buffer);

			if ((buffer_words.empty() || buffer_words.size() != 2)) eout("Invalid format. You should check " + file_name);
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



};

}


#endif /* CONTACT_STATE_READER_HPP */
