#ifndef FORCEFIELD_H_
#define FORCEFIELD_H_

#include "BaseObservable.h"

/**
 * @brief Prints ascii configurations. It can be extended to provide further
 * output types (e.g. for visualisation).
 *
 */

class ForceField: public BaseObservable {
protected:
	int only_type;
	std::set<int> visible_fields;
        std::set<int> hidden_fields;

	int Lx;
	int Ly;
	std::vector<number> f_field_x;
	std::vector<number> f_field_y;
	std::vector<number> f_field_passive_x;
	std::vector<number> f_field_passive_y;
	std::vector<number> f_field_active_x;
	std::vector<number> f_field_active_y;
	std::vector<number> phi_field;
	int delta = 8;
	int Lx_coarse;
	int Ly_coarse;
	std::vector<number> f_field_coarse_x;
	std::vector<number> f_field_coarse_y;
	std::vector<number> f_field_passive_coarse_x;
	std::vector<number> f_field_passive_coarse_y;
	std::vector<number> f_field_active_coarse_x;
	std::vector<number> f_field_active_coarse_y;
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
	virtual void calc_force_field(BaseField *p);
	//virtual void calc_stress_field(BaseField *p);

	/**
	 * @brief Returns the configuration output for the whole system. It does not comprise the headers.
	 *
	 * @param step
	 * @return
	 */
	virtual std::string force_field(llint step);

public:
	ForceField();
	virtual ~ForceField();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /* FORCEFIELD_H_ */

