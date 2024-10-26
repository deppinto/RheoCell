#ifndef ACTIVEMULTIFIELD_H_
#define ACTIVEMULTIFIELD_H_

#include "BaseInteraction.h"
#include "../Fields/MultiPhaseField.h"

/**
 * @brief Manages the interaction between simple patchy particles (as described in http://jcp.aip.org/resource/1/jcpsa6/v131/i1/p014504_s1)
 *
 * This interaction is selected with
 * interaction_type = simplemultifield
 *
 */

class ActiveMultiField: public BaseInteraction {
protected:

	number gamma;
	number lambda;
	number mu;
	int R;
	number kappa;
	number a0;
	number friction;
	number zetaS;
	number zetaS_active;

	/**
	 * @brief Patchy interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	number f_interaction(BaseField *p, int q);
	void calc_internal_forces(BaseField *p, int q);
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	void initFieldProperties(BaseField *p);
	std::vector<number> phi2;
	std::vector<number> sumS00;
	std::vector<number> sumS01;
	number velX;
	number velY;

public:
	ActiveMultiField();
	virtual ~ActiveMultiField();

	//void get_settings(input_file &inp) override;
	void init() override;
	void set_box(BaseBox *boxArg) override;

	void allocate_fields(std::vector<BaseField *> &fields) override;
	void apply_changes_after_equilibration() override;

	void begin_energy_computation() override;
	void begin_energy_computation(std::vector<BaseField *> &fields) override;
	void resetSums(int k) override;
	void updateFieldProperties(BaseField *p, int q, int k) override;

        void read_topology(std::vector<BaseField *> &fields) override;
	void check_input_sanity(std::vector<BaseField *> &fields) override;
	void get_settings(input_file &inp) override;
};

#endif /* ACTIVEMULTIFIELD_H_ */
