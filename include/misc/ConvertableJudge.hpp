#ifndef CONVERTABLE_JUDGE_HPP
#define CONVERTABLE_JUDGE_HPP

namespace cafemol {

namespace library {

bool isdigit(const std::string& String) {
	bool result = true;
	for (const char& Char : String) {
		result = result && (std::isdigit(static_cast<unsigned char>(Char)));
	}
	return result && (!String.empty());
}


bool isdigit(const std::vector<std::string>& Strings) {
	bool result = true;
	for (const std::string& String : Strings) {
		result = result && cafemol::library::isdigit(String);
	}
	return result && (!Strings.empty());
}

}
}

#endif /* CONVERTABLE_JUDGE_HPP */
