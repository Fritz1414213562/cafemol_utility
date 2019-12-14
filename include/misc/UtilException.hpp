#ifndef UTIL_EXCEPTION_HPP
#define UTIL_EXCEPTION_HPP
#include<string>

namespace cafemol {

namespace library {

class Util_Exception {

public:
	Util_Exception(const std::string& error_message) : what_exception(error_message) {}
	const std::string& what() {return what_exception;}

private:
	const std::string what_exception;

};
}
}

#endif /* UTIL_EXCEPTION_HPP */
