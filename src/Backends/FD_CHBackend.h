#ifndef FD_CHBACKEND_H_
#define FD_CHBACKEND_H_

#include "FDBackend.h"

/**
 * @brief Manages a FD simulation on CPU. It supports Predictor Corrector integration
 * This is the Cahn-Hilliard integration, i.e. it conserves the total field.
 */

class FD_CHBackend: public FDBackend {
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
	FD_CHBackend();
	virtual ~FD_CHBackend();

	void init();
	void get_settings(input_file &inp);
	void sim_step();
};

#endif /* FD_CHBACKEND_H_ */
