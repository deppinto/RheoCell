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
	//int R;
	number R1;
	number R2;
	number R_term;
	number Kg;
	number kappa;
	number a0;
	number friction;
	number friction_active;
	number zetaQ_self;
	number zetaQ_inter;
	number zetaQ_self_active;
	number zetaQ_inter_active;
	number J_Q;
	number J_Q_active;
	bool anchoring = false;
	number friction_cell;
	number friction_cell_active;
	number tolerance;
	number wall_slip; //if more than 1 there is no slip, if smaller than 1 there is slip. The length of slip on the walls is greater with decreasing values
	number passive_alpha;

	number F_total_x;
	number F_total_y;
	//std::vector<number> grad_free_energy_x;
	//std::vector<number> grad_free_energy_y;

	number f_interaction(BaseField *p, int q);
	void calc_internal_forces(BaseField *p, int q);
	void computeGlobalSums(BaseField *p, int q, bool update_global_sums=false);
	void initFieldProperties(BaseField *p);
	void update_anchoring(BaseField *p);

	Eigen::VectorXd vec_v_x;
	Eigen::VectorXd vec_f_x;
	//Eigen::SparseMatrix<double, Eigen::RowMajor> mat_m_x;
	Eigen::SparseMatrix<double> mat_m_x;
	Eigen::VectorXd vec_v_y;
	Eigen::VectorXd vec_f_y;
        void set_omp_tasks(int num_threads){Eigen::setNbThreads(num_threads);std::cout<<"TESTING: Set eigen openMP threads: "<<num_threads<<std::endl;};

	//Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  solverLU;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor> > solverCG;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solverCG;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::Lower|Eigen::Upper > solverCG;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solverCG;
	int restart_solver;
	int size_rows = 0;
	int size_rows_old = 0;
	int index, sub_q, other_site_patch, other_site_box;
	std::vector<int> neigh_values = std::vector<int> {0,1,2,3,5,6,7,8};
	//std::vector<number> weight_values = std::vector<number> {0.25, 0.5, 0.25, 0.5, 0., 0.5, 0.25, 0.5, 0.25};
	//std::vector<number> weight_values = std::vector<number> { 0.16667, 0.66667, 0.16667, 0.66667, 0., 0.66667, 0.16667, 0.66667, 0.16667};
	std::vector<number> weight_values = std::vector<number> {1, 1, 1, 1, 0, 1, 1, 1, 1};
	//std::vector<int> neigh_values = std::vector<int> {1,3,5,7};
	//std::vector<number> weight_values = std::vector<number> {0, 1, 0, 1, 0, 1, 0, 1, 0};
	//std::vector<int> neigh_values = std::vector<int> {4};
	std::vector<int> size_store_site_velocity_index;
	std::vector<int> store_site_velocity_index;
	std::vector<number> store_site_field;
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
	number get_total_force_x(BaseField *p, int q){return vec_v_x[q+field_start_index[p->index]] * friction;}
	number get_total_force_y(BaseField *p, int q){return vec_v_y[q+field_start_index[p->index]] * friction;}

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
