#ifndef BOXFACTORY_H_
#define BOXFACTORY_H_

#include "BaseBox.h"

/**
 * @brief Static factory class. Its only public method builds a {@link BaseBox}.
 *
 */
class BoxFactory {
public:
	BoxFactory() = delete;
	virtual ~BoxFactory() = delete;

	static BoxPtr make_box(input_file &inp);
};

#endif /* BOXFACTORY_H_ */
