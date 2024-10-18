#include "RCexception.h"
#include "Utils.h"

RCexception::RCexception(std::string ss, ...) {
	va_list ap;
	va_start(ap, ss);
	_error = Utils::sformat_ap(ss, ap);
	va_end(ap);
}

const char* RCexception::what() const throw() {
	return _error.c_str();
}
