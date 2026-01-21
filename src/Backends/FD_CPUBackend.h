#ifndef FD_CPUBACKEND_H_
#define FD_CPUBACKEND_H_

#include "FDBackend.h"

/**
 * @brief Manages a FD simulation on CPU. It supports Predictor Corrector integration.
 * This is the Allen-Cahn integration, i.e. does not conserve the total field.
 */

class FD_CPUBackend: public FDBackend {
protected:

	void first_step(bool store);
	void compute_forces();
	void second_step();

	number J0 = 1.0;
	bool analysis = false;

	void update_forces(BaseField *p, BaseField *q);
	vector<number> cos_x_table;
	vector<number> cos_y_table;
	vector<number> sin_x_table;
        vector<number> sin_y_table;

	number dphi, phi;
	std::vector<number> com_old = vector<number> {0., 0.};
	std::vector<number> temp = vector <number> {0., 0., 0., 0.};

public:
	FD_CPUBackend();
	virtual ~FD_CPUBackend();

	void init();
	void get_settings(input_file &inp);
	void sim_step();
};

#endif /* FD_CPUBACKEND_H_ */
