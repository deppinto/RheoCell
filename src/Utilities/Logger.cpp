#include "RCexception.h"
#include "Logger.h"

std::shared_ptr<Logger> Logger::logger = nullptr;

Logger::Logger() {
	log_open = false;
	log_stream = stderr;
	allow_log = true;
	debug = false;

	sprintf(log_level_strings[LOG_INFO], "%s", "INFO");
	sprintf(log_level_strings[LOG_WARNING], "%s", "WARNING");
	sprintf(log_level_strings[LOG_DEBUG], "%s", "DEBUG");
	sprintf(log_level_strings[LOG_ERROR], "%s", "ERROR");
}

Logger::~Logger() {
	if(log_open) {
		fclose(log_stream);
	}
}

void Logger::log(int log_level, const char *format, va_list &ap) {
	if(!allow_log)	{
		return;
	}

	if(log_level != LOG_NOTHING) {
		fprintf(log_stream, "%s: ", log_level_strings[log_level]);
	}
	vfprintf(log_stream, format, ap);
	va_end(ap);

	fprintf(log_stream, "\n");
	fflush(log_stream);
}

void Logger::set_stream(const char *filename) {
	FILE *buff = fopen(filename, "w");
	if(buff == nullptr) {
		throw RCexception("Log file '%s' is not writable", filename);
	}

	log_stream = buff;
	log_open = true;
}

void Logger::log(int log_level, const char *format, ...) {
	va_list ap;
	va_start(ap, format);
	log(log_level, format, ap);
}

void Logger::debug_log(const char *format, ...) {
	if(debug == true) {
		va_list ap;
		va_start(ap, format);
		log(LOG_DEBUG, format, ap);
	}
}

void Logger::get_settings(input_file &inp) {
	char filename[256];

	if(getInputString(&inp, "log_file", filename, 0) == KEY_FOUND) set_stream(filename);
	getInputBool(&inp, "debug", &debug, 0);
}

void Logger::init() {
	if(logger != nullptr) {
		throw RCexception("The logger has been already initialised");
	}

	logger = std::shared_ptr<Logger>(new Logger());
}

std::shared_ptr<Logger> Logger::instance() {
	if(logger == nullptr) {
		throw RCexception("Trying to access an uninitialised logger");
	}

	return logger;
}
