#ifndef WETMODEL_H_
#define WETMODEL_H_

#include "BaseInteraction.h"
#include "../Fields/MultiPhaseField.h"

/**
 * @brief Manages the interaction between fields with a dipolar force which aligns with shape
 *
 * This interaction is selected with
 * interaction_type = activenematic
 *
 */

class WetModel: public BaseInteraction {
protected:

	number gamma;
	number lambda;
	number omega;
	number mu;
	int R;
	number kappa;
	number a0;
	number friction;
	number zetaQ_self;
	number zetaQ_inter;
	number J_Q;
	bool anchoring = false;
	number friction_cell;

	number f_interaction(BaseField *p, int q);
	void calc_internal_forces(BaseField *p, int q);
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	void initFieldProperties(BaseField *p);
	void updateAnchoring(BaseField *p);
	std::vector<number> phi2;
	std::vector<number> sum_phi;
	std::vector<number> sumQ00;
	std::vector<number> sumQ01;
	std::vector<number> velocity_field_x;
	std::vector<number> velocity_field_y;
	number velX;
	number velY;

public:
	WetModel();
	virtual ~WetModel();

	//void get_settings(input_file &inp) override;
	void init() override;
	void set_box(BaseBox *boxArg) override;

	void allocate_fields(std::vector<BaseField *> &fields) override;

	void begin_energy_computation() override;
	void begin_energy_computation(std::vector<BaseField *> &fields) override;
	void resetSums(int k) override;
	void updateFieldProperties(BaseField *p, int q) override;

        void read_topology(std::vector<BaseField *> &fields) override;
	void check_input_sanity(std::vector<BaseField *> &fields) override;
	void get_settings(input_file &inp) override;
	void updateDirectedActiveForces(number dt, BaseField*p, bool store) override;
};

#endif /* WETMODEL_H_ */
