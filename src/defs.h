#ifndef DEFS_H_
#define DEFS_H_

#define PI 3.141592653589793238462643f
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define LRACOS(x) (((x) > 1) ? (number) 0 : ((x) < -1) ? (number) PI : acos(x))


#define CHECK_BOX(my_class, inp) 	std::string box_type ("");\
	if (getInputString(&inp, "box_type", box_type, 0) == KEY_FOUND) {\
		if (box_type.compare("Square") != 0) \
			throw ("%s only works with cubic box! Aborting", my_class);\
	}\


#define P_VIRTUAL (NULL)

#include "Utilities/LR_vector.h"
#include "Utilities/Logger.h"
#include "Utilities/parse_input/parse_input.h"

#include <string>
#include <memory>

using uint = uint32_t;
using llint = long long int;

#endif /* DEFS_H_ */
