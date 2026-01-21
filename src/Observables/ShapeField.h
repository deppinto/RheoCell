#ifndef SHAPEFIELD_H_
#define SHAPEFIELD_H_

#include "BaseObservable.h"

/**
 * @brief Prints ascii configurations. It can be extended to provide further
 * output types (e.g. for visualisation).
 *
 */

class ShapeField: public BaseObservable {
protected:
	int only_type;
	std::set<int> visible_fields;
        std::set<int> hidden_fields;

	int Lx;
	int Ly;
	std::vector<number> f_field_x;
	std::vector<number> f_field_y;
	std::vector<number> phi_field;
	int size_grid = 8;
	std::vector<number> f_field_coarse_x;
	std::vector<number> f_field_coarse_y;
	std::vector<number> phi_field_coarse;
	/**
	 * @brief Returns the configuration header(s)
	 *
	 * @param step
	 * @return
	 */
	virtual std::string headers(llint step);

	/**
	 * @brief Returns the portion of output for the given particle
	 *
	 * @param p
	 * @return
	 */
	virtual std::string field(BaseField *p);
	virtual void calc_field(BaseField *p);
	//virtual void calc_stress_field(BaseField *p);

	/**
	 * @brief Returns the configuration output for the whole system. It does not comprise the headers.
	 *
	 * @param step
	 * @return
	 */
	virtual std::string f_field(llint step);

public:
	ShapeField();
	virtual ~ShapeField();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /* SHAPEFIELD_H_ */

