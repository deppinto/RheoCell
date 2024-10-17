#ifndef HARDWALL_H_
#define HARDWALL_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines particles into a semispace. Given its "hard" nature, it may be used only iin Monte Carlo simulations.
 *
 */

class HardWall: public BaseForce {
public:
        HardWall();
        virtual ~HardWall();

	number position;
	number sigma;

	//std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual std::vector<number> value(llint step, std::vector<number> &pos);
	virtual number potential(llint step, std::vector<number> &pos);
};

#endif // HARDWALL_H
