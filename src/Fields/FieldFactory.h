#ifndef FIELDFACTORY_H_
#define FIELDFACTORY_H_

#include "BaseField.h"

/**
 * @brief Static factory class. Its only public method builds a {@link BaseBox}.
 *
 * @verbatim
 [list_type = multiphasefield (Type of simulation box for CPU simulations.)]
 @endverbatim
 */
class FieldFactory {
public:
	FieldFactory() = delete;
	virtual ~FieldFactory() = delete;

	static FieldPtr make_field();
};

#endif /* FIELDFACTORY_H_ */
