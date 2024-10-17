#ifndef FD_CPUBACKEND_H_
#define FD_CPUBACKEND_H_

#include "FDBackend.h"

/**
 * @brief Manages a FD simulation on CPU. It supports Predictor Corrector integration
 */

class FD_CPUBackend: public FDBackend {
protected:

	void first_step(bool store);
	void first_eq_step(bool store);
	void compute_forces();
	void second_step();

	number J0= 1.0;

	void update_forces(BaseField *p, BaseField *q);
	vector<number> cos_x_table;
	vector<number> cos_y_table;
	vector<number> sin_x_table;
        vector<number> sin_y_table;


public:
	FD_CPUBackend();
	virtual ~FD_CPUBackend();

	void init();
	void get_settings(input_file &inp);
	void sim_step();
	void eq_step();
};

#endif /* FD_CPUBACKEND_H_ */
