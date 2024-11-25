#ifndef SIMPLEMULTIFIELD_H_
#define SIMPLEMULTIFIELD_H_

#include "BaseInteraction.h"
#include "../Fields/MultiPhaseField.h"

/**
 * @brief Manages the interaction between simple patchy particles (as described in http://jcp.aip.org/resource/1/jcpsa6/v131/i1/p014504_s1)
 *
 * This interaction is selected with
 * interaction_type = simplemultifield
 *
 */

class SimpleMultiField: public BaseInteraction {
protected:

	number gamma;
	number lambda;
	number mu;
	//int R;
	number kappa;
	number a0;

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
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	std::vector<number> phi2;

public:
	SimpleMultiField();
	virtual ~SimpleMultiField();

	//void get_settings(input_file &inp) override;
	void init() override;
	void set_box(BaseBox *boxArg) override;

	void allocate_fields(std::vector<BaseField *> &fields) override;

	void begin_energy_computation() override;
	void begin_energy_computation(std::vector<BaseField *> &fields) override;
	void resetSums(int k) override;

        void read_topology(std::vector<BaseField *> &fields) override;
	void check_input_sanity(std::vector<BaseField *> &fields) override;
	void get_settings(input_file &inp) override;
};

#endif /* SIMPLEMULTIFIELD_H_ */
