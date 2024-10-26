#ifndef FDBACKEND_H_
#define FDBACKEND_H_

#include "SimBackend.h"

using namespace std;

/**
 * @brief Abstract class for Predictor Corrector simulations.
 *
 * This class sets up a basic FD simulation. It does not do any sensible computation but
 * takes care of the most basic input/output operations associated with MD simulations.
 *
 */

class FDBackend: public SimBackend {
protected:
	number dt = 0.;

	// timers
	std::shared_ptr<Timer> timer_first_step;
	std::shared_ptr<Timer> timer_forces;
	std::shared_ptr<Timer> timer_testing;

	bool lees_edwards;
	number shear_rate;

public:
	FDBackend();
	virtual ~FDBackend();

	void get_settings(input_file &inp);
	void init();

	virtual void print_observables();
};

#endif /* MDBACKEND_H_ */
