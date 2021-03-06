#ifndef CLUSTERING_RESULT_READER_HPP
#define CLUSTERING_RESULT_READER_HPP
#include<OtherFormat/ClusteringResultParser.hpp>
#include<iostream>
#include<vector>
#include<array>
#include<string>
#include<fstream>


namespace cafemol {

class ClusteringResultReader : public ClusteringResultParser {

public:

	ClusteringResultReader() : ClusteringResultParser() {}
	ClusteringResultReader(const std::string& input_filename) : ClusteringResultParser(input_filename) {}

	std::vector<std::vector<std::size_t>> read_ClassifiedIDs() {

		if (!input_file.is_open()) eout("No file has loaded.");

		std::string buffer;
		std::vector<std::vector<std::size_t>> result;
		bool is_on_cluster_ids_block = false;

		while (std::getline(input_file, buffer)) {

			if (buffer.empty()) continue;

			std::vector<std::string> buffer_words = split_Line(buffer);

			if ((!is_on_cluster_ids_block) && (buffer_words.size() < 2)) continue;

			if ((is_on_cluster_ids_block) && (buffer_words[0] == ">>>>")) {
				is_on_cluster_ids_block = false;
				continue;
			}

			if ((buffer_words[0] == "<<<<") && (buffer_words[1] == "cluster_ids")) {
				is_on_cluster_ids_block = true;
				continue;
			}
			if (!is_on_cluster_ids_block) continue;

			if (buffer_words[0] == "Cluster:") continue;
			result.push_back(split_ClusterIDLine(buffer));
		}

		close_File();
		open_File();

		return result;
	}

	float read_SimilaritySum() {

		if (!input_file.is_open()) eout("No file has loaded.");

		std::string buffer;
		float result;
		bool is_SS_read = false;
		bool is_on_SS_block = false;

		while (std::getline(input_file, buffer)) {
			if (buffer.empty()) continue;

			std::vector<std::string> buffer_words = split_Line(buffer);

			if ((!is_on_SS_block) && (buffer_words.size() < 2)) continue;

			if ((is_on_SS_block) && (buffer_words[0] == ">>>>")) {
				is_on_SS_block = false;
				continue;
			}

			if ((buffer_words[0] == "<<<<") && (buffer_words[1] == "similarity_sum")) {
				is_on_SS_block = true;
				continue;
			}
			if (!is_on_SS_block) continue;

			if ((buffer_words[0] == "SS") && (buffer_words[1] == "=")) {
				result = std::stof(buffer_words[2]);
				is_SS_read = true;
			}
			else eout("Invalid format on 'similarity_sum' block");
		}
		if (!is_SS_read) eout("This result file might be broken because 'similarity_sum' block could not be found");

		close_File();
		open_File();

		return result;
	}

};
}


#endif /* CLUSTERING_RESULT_READER_HPP */
