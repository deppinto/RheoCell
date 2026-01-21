#ifndef DIFFERENTIALADHESION_H_
#define DIFFERENTIALADHESION_H_

#include "BaseInteraction.h"
#include "../Fields/MultiPhaseField.h"

/**
 * @brief Manages the interaction between simple patchy particles (as described in http://jcp.aip.org/resource/1/jcpsa6/v131/i1/p014504_s1)
 *
 * This interaction is selected with
 * interaction_type = simplemultifield
 *
 */

class DifferentialAdhesion: public BaseInteraction {
protected:

	number gamma;
	number lambda;
	number mu;
	//int R;
	number kappa;
	number omega;
	number a0;
	number zetaQ_self;
	number zetaQ_inter;
	number zetaQ_self_active;
	number zetaQ_inter_active;
	number strain_rate;
	number strain_rate_active;

	int size_rows = 0;
	int size_rows_old = 0;
	std::vector<number> vec_omega;
	std::vector<number> phi_omega;
	int store_max_size = 0;
	std::vector<int> size_store_site_omega_index;
	std::vector<int> store_site_omega_index;
	std::vector<int> store_site_omega_sub;
	std::vector<number> store_site_field;
	std::vector<int> field_start_index;
	/**
	 * @brief Patchy interaction between two particles.
	 *
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	//number f_interaction(BaseField *p, int q);
	number f_interaction(BaseField *p, int q, BaseField *pp, int qq);
	void construct_omega(std::vector<BaseField *> &fields);
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	void calc_internal_forces(BaseField *p, int q);
	void initFieldProperties(BaseField *p);
	std::vector<number> phi2;
	std::vector<number> sumQ00;
	std::vector<number> sumQ01;

public:
	DifferentialAdhesion();
	virtual ~DifferentialAdhesion();

	//void get_settings(input_file &inp) override;
	void init() override;
	void set_box(BaseBox *boxArg) override;

	void allocate_fields(std::vector<BaseField *> &fields) override;

	void begin_energy_computation() override;
	void begin_energy_computation(std::vector<BaseField *> &fields) override;
	void resetSums(int k) override;
	void updateFieldProperties(BaseField *p, int q, int k) override;

        void read_topology(std::vector<BaseField *> &fields) override;
	void check_input_sanity(std::vector<BaseField *> &fields) override;
	void get_settings(input_file &inp) override;
	void apply_changes_after_equilibration() override;
	void updateDirectedActiveForces(number dt, BaseField*p, bool store) override;
};

#endif /* DIFFERENTIALADHESION_H_ */
