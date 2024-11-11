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
#include <omp.h>
#include <thread>

using uint = uint32_t;
using llint = long long int;

/*!
omp Template: loop over the function with omp or not
the syntax requires that the first argument of the function is the "index" of whatever the function acts on.
So, if the function is f(int, double, double,Index2D,...) then this template function should be called by:
ompFunctionLoop(ompThreadNum,maxIdx, f, double, double,Index2D,...).
*/
template< typename... Args>
void ompFunctionLoop(int nThreads, int maxIdx, void (*fPointer)(int, Args...), Args... args)
    {
    if(nThreads <= 1)
        {
        for(int idx = 0; idx < maxIdx; ++idx)
            fPointer(idx,args...);
        }
    else
        {
	    #pragma omp parallel for num_threads(nThreads)
        for(int idx = 0; idx < maxIdx; ++idx)
            fPointer(idx,args...);
        }
    };


#endif /* DEFS_H_ */
