#ifndef LOGGER_H_
#define LOGGER_H_

#define OX_LOG Logger::instance()->log
#define OX_DEBUG Logger::instance()->debug_log

#include <cstdarg>
#include <memory>

#include "parse_input/parse_input.h"

/**
 * @brief Simple logger
 *
 * This logger is a singleton (http://www.yolinux.com/TUTORIALS/C++Singleton.html) used to log
 * infos, warnings and errors either on a log file or to stderr.
 */
class Logger {
private:
	/// Static pointer to the only allowed Logger instance
	static std::shared_ptr<Logger> logger;
	/// If true, the logger prints also debug informations
	bool debug;
	/// If false, the logger won't print anything
	bool allow_log;

	/// Output stream
	FILE *log_stream;
	/// Output string for each log level
	char log_level_strings[4][256];
	/// If true, the log stream can be used
	bool log_open;

	void log(int log_level, const char *format, va_list &ap);
	void set_stream(const char *filename);

	/**
	 * @brief Default constructor. It is kept private to enforce the singleton pattern.
	 *
	 * @return
	 */
	Logger();

public:
	enum log_levels {
		LOG_INFO = 0,
		LOG_WARNING = 1,
		LOG_DEBUG = 2,
		LOG_ERROR = 3,
		LOG_NOTHING = 4
	};

	Logger(Logger const&) = delete;

	/**
	 * @brief Read options from the input file
	 *
	 * @param inp reference to an input_file struct containing a dictionary of key/values
	 */
	void get_settings(input_file &inp);

	/**
	 * @brief Initialize the logger. It is static to avoid more than one instantiation
	 */
	static void init();

	/**
	 * @brief Variadic method. Does the actual logging
	 *
	 * @param log_level it should be one of the log_levels
	 * @param format printf-like string
	 */
	void log(int log_level, const char *format, ...);

	/**
	 * @brief Variadic method. Used to log debug messages.
	 *
	 * @param format
	 */
	void debug_log(const char *format, ...);

	/**
	 * @brief Disables logging.
	 */
	void disable_log() { allow_log = false; }

	/**
	 * @brief Enabled logging.
	 */
	void enable_log() { allow_log = true; }

	/**
	 * @brief Returns the pointer to the log stream used by this Logger
	 *
	 * @return pointer to the log stream (it can be stdout, stderr or an open file)
	 */
	FILE *get_log_stream() { return log_stream; }

	/**
	 * @brief Returns the actual logger. Static method to enforce the singleton pattern.
	 *
	 * @return Pointer to an already initialized logger
	 */
	static std::shared_ptr<Logger> instance();

	virtual ~Logger();
};

#endif /* LOGGER_H_ */
