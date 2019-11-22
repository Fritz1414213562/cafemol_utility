#include"../../include/misc/ErrorMessage.hpp"

int main() {

	cafemol::error_handling::Error_Output eo = cafemol::error_handling::Error_Output();
	eo("This is a test.", "But this is not just a test.", "Do you think so?");

	return 0;
}
