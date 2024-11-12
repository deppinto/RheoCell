#ifndef WETMODEL_H_
#define WETMODEL_H_

#include "BaseInteraction.h"
#include "../Fields/MultiPhaseField.h"
#include <Eigen/Sparse>
#include <Eigen/Core>

/**
 * @brief Manages the interaction between fields with a dipolar force which aligns with shape and introduces viscosity
 *
 * This interaction is selected with
 * interaction_type = wetmodel
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
	number friction_active;
	number zetaQ_self;
	number zetaQ_inter;
	number zetaQ_self_active;
	number zetaQ_inter_active;
	number J_Q;
	bool anchoring = false;
	number friction_cell;

	number f_interaction(BaseField *p, int q);
	void calc_internal_forces(BaseField *p, int q);
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	void initFieldProperties(BaseField *p);
	void update_anchoring(BaseField *p);

	Eigen::VectorXd vec_v_x;
	Eigen::VectorXd vec_f_x;
	Eigen::SparseMatrix<double, Eigen::RowMajor> mat_m_x;
	//Eigen::SparseMatrix<double> mat_m_x;
	Eigen::VectorXd vec_v_y;
	Eigen::VectorXd vec_f_y;
        void set_omp_tasks(int num_threads){Eigen::setNbThreads(num_threads);std::cout<<"TESTING: Set eigen openMP threads: "<<num_threads<<std::endl;};

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  solverLU;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor> > solverCG;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solverCG;
	int size_rows = 0;
	int index, sub_q, other_site_patch, other_site_box;
	std::vector<int> neigh_values = std::vector<int> {5,3,1,7};
	std::vector<int> size_store_site_velocity_index;
	std::vector<int> store_site_velocity_index;
	std::vector<int> field_start_index;
	int store_max_size = 0;
	std::vector<number> phi2;
	std::vector<number> sumQ00;
	std::vector<number> sumQ01;
	std::vector<number> sum_phi;
	number velX;
	number velY;

public:
	WetModel();
	virtual ~WetModel();

	//void get_settings(input_file &inp) override;
	void init() override;
	void set_box(BaseBox *boxArg) override;

	void allocate_fields(std::vector<BaseField *> &fields) override;
	void apply_changes_after_equilibration() override;
	number get_velocity_x(BaseField *p, int q){return vec_v_x[q+field_start_index[p->index]];}
	number get_velocity_y(BaseField *p, int q){return vec_v_y[q+field_start_index[p->index]];}
	//number get_velocity_x(BaseField *p, int q){return vec_v_x[q+p->index*p->subSize];}
	//number get_velocity_y(BaseField *p, int q){return vec_v_y[q+p->index*p->subSize];}

	void begin_energy_computation() override;
	void begin_energy_computation(std::vector<BaseField *> &fields) override;
	void resetSums(int k) override;
	void updateFieldProperties(BaseField *p, int q, int k) override;

        void read_topology(std::vector<BaseField *> &fields) override;
	void check_input_sanity(std::vector<BaseField *> &fields) override;
	void get_settings(input_file &inp) override;
	void updateDirectedActiveForces(number dt, BaseField*p, bool store) override;
};

#endif /* WETMODEL_H_ */
