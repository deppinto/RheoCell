#ifndef LEBCACTIVENEMATIC_H_
#define LEBCACTIVENEMATIC_H_

#include "BaseInteraction.h"
#include "../Fields/MultiPhaseField.h"

/**
 * @brief Manages the interaction between fields with a dipolar force which aligns with shape
 *
 * This interaction is selected with
 * interaction_type = wetmodel
 *
 */

class LEBcActiveNematic: public BaseInteraction {
protected:

	number gamma;
	number lambda;
	number omega;
	number mu;
	number kappa;
	number a0;
	number friction;
	number zetaQ_self;
	number zetaQ_inter;
	number zetaQ_self_active;
	number zetaQ_inter_active;
	number J_Q;
	number shear_rate;
	bool anchoring = false;

	number f_interaction(BaseField *p, int q);
	void calc_internal_forces(BaseField *p, int q);
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	void initFieldProperties(BaseField *p);
	void updateAnchoring(BaseField *p);
	std::vector<number> phi2;
	std::vector<number> sumQ00;
	std::vector<number> sumQ01;
	number velX;
	number velY;

public:
	LEBcActiveNematic();
	virtual ~LEBcActiveNematic();

	//void get_settings(input_file &inp) override;
	void init() override;
	void set_box(BaseBox *boxArg) override;

	void allocate_fields(std::vector<BaseField *> &fields) override;
	void apply_changes_after_equilibration() override;

	void begin_energy_computation() override;
	void begin_energy_computation(std::vector<BaseField *> &fields) override;
	void resetSums(int k) override;
	void updateFieldProperties(BaseField *p, int q, int k) override;
	number get_velocity_x(BaseField *p, int q){
		int y = p->map_sub_to_box[q]/box->getXsize();
		if(y==0) return p->velocityX[q] - shear_rate * box->getYsize();
		if(y==box->getYsize()-1) return p->velocityX[q] + shear_rate * box->getYsize();
		else return p->velocityX[q];
	}
	number get_velocity_y(BaseField *p, int q){return p->velocityY[q];}

        void read_topology(std::vector<BaseField *> &fields) override;
	void check_input_sanity(std::vector<BaseField *> &fields) override;
	void get_settings(input_file &inp) override;
	void updateDirectedActiveForces(number dt, BaseField*p, bool store) override;
};

#endif /* LEBCACTIVENEMATIC_H_ */
