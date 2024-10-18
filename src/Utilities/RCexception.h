#ifndef RCEXCEPTION_H_
#define RCEXCEPTION_H_

#include <string>
#include <cstdarg>

/**
 * @brief Generic oxDNA exception. It should be handled by either SimManager or AnalysisManager.
 */
class RCexception : public std::exception {
private:
	std::string _error;
public:
	RCexception(std::string ss, ...);
	// I'm not sure why, but the throw() bit is needed to avoid a 'looser throw specifier' error
	virtual ~RCexception() throw() {};

	/**
	 * @brief Returns the actual error message.
	 *
	 * @return error message associated to the exception
	 */
	virtual const char* what() const throw();
};

#endif /* RCEXCEPTION_H_ */
