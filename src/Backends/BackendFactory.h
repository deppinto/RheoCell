#ifndef BACKENDFACTORY_H_
#define BACKENDFACTORY_H_

#include "SimBackend.h"

/**
 * @brief Static factory class. It exposes a single static method that builds a {@link SimBackend simulation backend}.
 *
 * This class can not be instantiated. It provides a single method that
 * parses the input file and builds the right simulation backend.
 *
 */

class BackendFactory {
public:
	BackendFactory() = delete;
	virtual ~BackendFactory() = delete;

	/**
	 * @brief Builds the backend.
	 *
	 * @param inp
	 * @return a pointer to the newly created backend
	 */
	static std::shared_ptr<SimBackend> make_backend(input_file &inp);
};

#endif /* BACKENDFACTORY_H_ */
