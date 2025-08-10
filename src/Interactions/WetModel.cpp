#include "WetModel.h"
#include "../Utilities/Utils.h"

WetModel::WetModel() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				omega(0.004),
				mu(3.),
				Kg(0.1),
				kappa(0.1),
				friction(1.),
				zetaQ_self(0),
				zetaQ_inter(0),
				J_Q(0),
				friction_cell(0.),
				tolerance(0.0001),
				wall_slip(0.5),
				passive_alpha(1.){
	a0=PI*R*R;
}

WetModel::~WetModel() {

}


void WetModel::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "R", &R, 0);
	getInputNumber(&inp, "lambda", &lambda, 0);
	getInputNumber(&inp, "omega", &omega, 0);
	getInputNumber(&inp, "gamma", &gamma, 0);
	getInputNumber(&inp, "mu", &mu, 0);
	getInputNumber(&inp, "kappa", &kappa, 0);
	getInputNumber(&inp, "zetaQ_self", &zetaQ_self_active, 0);
	getInputNumber(&inp, "zetaQ_inter", &zetaQ_inter_active, 0);
	getInputNumber(&inp, "J_Q", &J_Q_active, 0);
	getInputBool(&inp, "anchoring", &anchoring, 0);
	getInputNumber(&inp, "friction_cell", &friction_cell_active, 0);
	getInputNumber(&inp, "friction", &friction_active, 0);
	getInputNumber(&inp, "CGtolerance", &tolerance, 0);
	getInputNumber(&inp, "wall_slip", &wall_slip, 0);
	getInputNumber(&inp, "Kg", &Kg, 0);
	getInputNumber(&inp, "passive_alpha", &passive_alpha, 0);
}

void WetModel::init() {
        a0=PI*R*R;
	/*number R_eff = R;
	number AR = 3.;
	R1 = (R_eff * sqrt(AR)) * (R_eff * sqrt(AR));
	R2 = (R_eff / sqrt(AR)) * (R_eff / sqrt(AR));
	R_term = (1/R1) - (1/R2);*/
	//R1 = 1;
	//R2 = 1;
	store_max_size=20;
	if(tolerance>0)
		solverCG.setTolerance(tolerance);
        set_omp_tasks(omp_thread_num);
	restart_solver=0;
	std::cout<<"TESTING: Running with tolerance (CG): "<<tolerance<<std::endl;
}

void WetModel::read_topology(std::vector<BaseField*> &fields) {
        int N = fields.size();
	field_start_index.resize(N);

        std::ifstream topology(topology_filename, std::ios::in);
        if(!topology.good()) {
                throw RCexception("Can't read topology file '%s'. Aborting", topology_filename);
        }

        allocate_fields(fields);
        for(int i = 0; i < N; i++) {
                fields[i]->index = i;
		fields[i]->get_interaction_values(R);
        }
}

void WetModel::allocate_fields(std::vector<BaseField *> &fields) {
        for(int i = 0; i < (int) fields.size(); i++) {
                fields[i] = new MultiPhaseField();
        }
}

void WetModel::apply_changes_after_equilibration(){

	BaseInteraction::apply_changes_after_equilibration();
	zetaQ_self=zetaQ_self_active;
	zetaQ_inter=zetaQ_inter_active;
	J_Q=J_Q_active;
	friction=friction_active;
	friction_cell=friction_cell_active;
}

void WetModel::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
	if(box->lees_edwards)throw RCexception("Interaction is not compatible with LEBc. Aborting");
	phi2.resize(Lx*Ly);
	sumQ00.resize(Lx*Ly);
	sumQ01.resize(Lx*Ly);
	sum_phi.resize(Lx*Ly);
	//grad_free_energy_x.resize(Lx*Ly);
	//grad_free_energy_y.resize(Lx*Ly);
	size_store_site_velocity_index.resize(Lx*Ly);
	store_site_velocity_index.resize(Lx*Ly*store_max_size);
	store_site_field.resize(Lx*Ly*store_max_size);
	for(int i =0; i<Lx*Ly; i++){resetSums(i);size_store_site_velocity_index[i]=0;}
}

void WetModel::resetSums(int k) {
	phi2[k]=0;
        sumQ00[k]=0;
        sumQ01[k]=0;
        sum_phi[k]=0;
	//grad_free_energy_x[k]=0.;
	//grad_free_energy_y[k]=0.;
}

void WetModel::initFieldProperties(BaseField *p) {

	//number ytop, ybottom, xleft, xright;
	for(int q=0; q<p->subSize;q++) {
		int k = p->GetSubIndex(q, box);
		BaseInteraction::updateFieldProperties(p, q, k);

		/*if(p->neighbors_sub[5+q*9]==-1)xright=p->fieldScalar[q];
		else xright=p->fieldScalar[p->neighbors_sub[5+q*9]];
		if(p->neighbors_sub[3+q*9]==-1)xleft=p->fieldScalar[q];
		else xleft=p->fieldScalar[p->neighbors_sub[3+q*9]];
		if(p->neighbors_sub[7+q*9]==-1)ytop=p->fieldScalar[q];
		else ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
		if(p->neighbors_sub[1+q*9]==-1)ybottom=p->fieldScalar[q];
		else ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
		if(p->neighbors_sub[5+q*9]==-1)xright=0;
		else xright=p->fieldScalar[p->neighbors_sub[5+q*9]];
		if(p->neighbors_sub[3+q*9]==-1)xleft=0;
		else xleft=p->fieldScalar[p->neighbors_sub[3+q*9]];
		if(p->neighbors_sub[7+q*9]==-1)ytop=0;
		else ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
		if(p->neighbors_sub[1+q*9]==-1)ybottom=0;
		else ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
	        number dx = .5*( xright - xleft );
	        number dy = .5*( ytop - ybottom );*/

	        //number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
	        //number dy = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
	        number dx = BaseInteraction::derivX(p, q, k);
	        number dy = BaseInteraction::derivY(p, q, k);
	        p->fieldDX[q] = dx;
	        p->fieldDY[q] = dy;

		p->S00 += -0.5*(dx*dx-dy*dy);
		p->S01 += -dx*dy;
		p->NQ00 = 0.;
		p->NQ01 = 0.;
	}
}


void WetModel::updateFieldProperties(BaseField *p, int q, int k) {
	/*BaseInteraction::updateFieldProperties(p, q, k);
	number dx = p->fieldDX[q]; 
	number dy = p->fieldDY[q]; 
	p->S00 += -0.5*(dx*dx-dy*dy);
	p->S01 += -dx*dy;*/

	p->S00 += -0.5*(p->fieldDX[q]*p->fieldDX[q]-p->fieldDY[q]*p->fieldDY[q]);
	p->S01 += -p->fieldDX[q]*p->fieldDY[q];
	p->NQ00 = 0.;
	p->NQ01 = 0.;
}


void WetModel::check_input_sanity(std::vector<BaseField *> &fields) {

}


void WetModel::begin_energy_computation() {
		
        for(int i = 0; i < CONFIG_INFO->N(); i++) {
                initFieldProperties(CONFIG_INFO->fields()[i]);
        }
}

void WetModel::begin_energy_computation(std::vector<BaseField *> &fields) {

	std::fill(size_store_site_velocity_index.begin(), size_store_site_velocity_index.end(), 0);
	if(restart_solver<10)
		size_rows_old = size_rows;
	else{
		size_rows_old = 0;
		restart_solver = 0;
	}

	size_rows=0;
	for(auto p : fields) {
		field_start_index[p->index]=size_rows;
		size_rows += p->subSize;
		for(int q=0; q<p->subSize;q++)
			computeGlobalSums(p, q, false);
	}
	if(size_rows != size_rows_old){
		vec_v_x.resize(size_rows);
		vec_f_x.resize(size_rows);
		vec_v_y.resize(size_rows);
		vec_f_y.resize(size_rows);
		mat_m_x.resize(size_rows, size_rows);
		mat_m_x.setZero();
	}
	else
		mat_m_x.setZero();


	//F_total_x = 0.;
	//F_total_y = 0.;
	//number grad_x = 0.;
	//number grad_y = 0.;
        U = (number) 0;
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++)
			U += f_interaction(p, q);

                /*for(int q=0; q<p->subSize;q++){
        		grad_free_energy_x[p->map_sub_to_box[q]] += (-1) * 0.5 * p->fieldScalar[q] * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] );
        		grad_free_energy_y[p->map_sub_to_box[q]] += (-1) * 0.5 * p->fieldScalar[q] * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] );
        		grad_x += (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] );
        		grad_y += (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] );
		}*/
        }
	//for(int q=0; q<box->getXsize()*box->getYsize();q++){grad_x+=grad_free_energy_x[q];grad_y+=grad_free_energy_y[q];}
	//std::cout<<"Passive Forces: "<<F_total_x<<" "<<F_total_y<<" "<<grad_x<<" "<<grad_y<<" "<<grad_free_energy_x[0]<<" "<<grad_free_energy_y[0] <<std::endl;


        K = (number) 0;
	std::vector<Eigen::Triplet<double>> tri_t_x;
	//F_total_x = 0.;
	//F_total_y = 0.;
	//mat_m_x.setZero();
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++){
			calc_internal_forces(p, q);

			if(friction_cell==0){
				vec_v_x[q+field_start_index[p->index]] = vec_f_x[q+field_start_index[p->index]]/friction;
				vec_v_y[q+field_start_index[p->index]] = vec_f_y[q+field_start_index[p->index]]/friction;
				continue;
			}

			//populate sparse matrix
			if(box->getWalls(p->map_sub_to_box[q])<1.)//wall_slip)
				//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction+4*friction_cell)));
				//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction+friction_cell)));
				//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction+8*friction_cell)));
				tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction)));
			else
				tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], 1.0));
			//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction) ));

			/*other_site_patch = q;
			other_site_box = p->map_sub_to_box[other_site_patch];
			if(box->getWalls(p->map_sub_to_box[q])<wall_slip)
				tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], other_site_patch+field_start_index[p->index] , (double)(-4*friction_cell)));
			if(sum_phi[other_site_box]!=0){
				for(int i=0; i<size_store_site_velocity_index[other_site_box]; i++){
					if(box->getWalls(p->map_sub_to_box[q])<wall_slip)
						tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(4*friction_cell*p->fieldScalar[other_site_patch]/sum_phi[other_site_box])));
				}
			}*/


			for(auto j : neigh_values){
				//other_site_patch = p->neighbors_sub[j+q*9];
				//other_site_box = p->map_sub_to_box[other_site_patch];
				other_site_box = box->neighbors[j+p->map_sub_to_box[q]*9];
				//if(sum_phi[p->map_sub_to_box[q]]==0)continue;
				//if(sum_phi[other_site_box]==0)continue;

				//if(box->getWalls(p->map_sub_to_box[q])<wall_slip)
					//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], other_site_patch+field_start_index[p->index] , (double)(friction_cell)));

				for(int i=0; i<size_store_site_velocity_index[other_site_box]; i++){
					//number weight = 1.0; //p->fieldScalar[q] * store_site_field[i+other_site_box*store_max_size]/(sum_phi[p->map_sub_to_box[q]] * sum_phi[other_site_box]);
					//tri_t_x.push_back(Eigen::Triplet<double> (store_site_velocity_index[i+other_site_box*store_max_size], q+field_start_index[p->index], (double)(-friction_cell*p->fieldScalar[other_site_patch]/sum_phi[other_site_box])));
					
					//if(store_site_velocity_index[i+other_site_box*store_max_size]>=field_start_index[p->index] && store_site_velocity_index[i+other_site_box*store_max_size]<field_start_index[p->index]+p->subSize)continue;

					if(box->getWalls(p->map_sub_to_box[q])<wall_slip){
						tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction_cell * weight_values[j])));
						tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell * weight_values[j])));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction_cell)));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell)));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction_cell)*weight));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell)*weight));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell*store_site_field[i+other_site_box*store_max_size]/sum_phi[other_site_box])));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell*p->fieldScalar[q]/sum_phi[p->map_sub_to_box[q]])));
					}

					//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(friction_cell*p->fieldScalar[other_site_patch]/sum_phi[other_site_box])));
				}
			}
		}
        }
	//std::cout<<"Active Forces: "<<F_total_x<<" "<<F_total_y<<std::endl;

	/*
	//std::cout<<"start solver: "<<size_rows/2 <<std::endl;
	mat_m_x.setFromTriplets(tri_t_x.begin(), tri_t_x.end());
	mat_m_x.makeCompressed();
	//std::cout<<"copy done"<<std::endl;
	solverLU.analyzePattern(mat_m_x);
	solverLU.factorize(mat_m_x);
	vec_v_x = solverLU.solve(vec_f_x);
	vec_v_y = solverLU.solve(vec_f_y);
	*/

	if(friction_cell!=0){
	//std::cout<<Eigen::nbThreads()<<std::endl;
	//std::cout<<"start solver: "<<size_rows/2 <<std::endl;
		mat_m_x.setFromTriplets(tri_t_x.begin(), tri_t_x.end());
		solverCG.compute(mat_m_x);
		if(size_rows != size_rows_old){
			vec_v_x = solverCG.solve(vec_f_x);
			vec_v_y = solverCG.solve(vec_f_y);
			restart_solver = 0;
		}
		else{
			vec_v_x = solverCG.solveWithGuess(vec_f_x, vec_v_x);
			vec_v_y = solverCG.solveWithGuess(vec_f_y, vec_v_y);
			restart_solver+=1;
		}
	}
	//std::cout << "#iterations:     " << solverCG.iterations() << std::endl;
	//std::cout << "estimated error: " << solverCG.error()      << std::endl;	
	//std::cout<<"end solver"<<std::endl;

	/*number V_All_x=0.;
	number V_All_y=0.;
	number V_CoM_x=0.;
	number V_CoM_y=0.;
	number V_gammaphi_x=0.;
	number V_gammaphi_y=0.;
	number F_CoM_x=0.;
	number F_CoM_y=0.;
	number V_total_x = 0.;
	number V_total_y = 0.;
	number V_phi_x=0.;
	number V_phi_y=0.;
	number F_phi_x=0.;
	number F_phi_y=0.;
	F_total_x = 0.;
	F_total_y = 0.;
	//for(int z=0; z<size_rows;z++){
	//	F_total_x += vec_v_x[z];
	//	F_total_y += vec_v_y[z];
	//}
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++){
			//F_total_x += friction * vec_v_x[q+field_start_index[p->index]] + 8 * vec_v_x[q+field_start_index[p->index]];
			//F_total_y += friction * vec_v_y[q+field_start_index[p->index]] + 8 * vec_v_y[q+field_start_index[p->index]];
			V_total_x += friction * vec_v_x[q+field_start_index[p->index]];
			V_total_y += friction * vec_v_y[q+field_start_index[p->index]];
			V_All_x += vec_v_x[q+field_start_index[p->index]];
			V_All_y += vec_v_y[q+field_start_index[p->index]];
			V_gammaphi_x += vec_v_x[q+field_start_index[p->index]] * p->fieldDX[q];
			V_gammaphi_y += vec_v_y[q+field_start_index[p->index]] * p->fieldDY[q];
			V_phi_x += p->fieldScalar[q] * vec_v_x[q+field_start_index[p->index]];
			V_phi_y += p->fieldScalar[q] * vec_v_y[q+field_start_index[p->index]];

			F_total_x += vec_f_x[q+field_start_index[p->index]];
			F_total_y += vec_f_y[q+field_start_index[p->index]];
			F_phi_x += p->fieldScalar[q] * vec_f_x[q+field_start_index[p->index]];
			F_phi_y += p->fieldScalar[q] * vec_f_y[q+field_start_index[p->index]];

			V_CoM_x += p->fieldScalar[q] * vec_v_x[q+field_start_index[p->index]] / p->area;
			V_CoM_y += p->fieldScalar[q] * vec_v_y[q+field_start_index[p->index]] / p->area;
			F_CoM_x += p->fieldScalar[q] * vec_f_x[q+field_start_index[p->index]] / p->area;
			F_CoM_y += p->fieldScalar[q] * vec_f_y[q+field_start_index[p->index]] / p->area;
			for(auto j : neigh_values){
				other_site_box = box->neighbors[j+p->map_sub_to_box[q]*9];
				for(int i=0; i<size_store_site_velocity_index[other_site_box]; i++){
					//F_total_x += (double)(-friction_cell*p->fieldScalar[q]/sum_phi[p->map_sub_to_box[q]]) * vec_v_x[q+store_site_velocity_index[i+other_site_box*store_max_size]];
					//F_total_y += (double)(-friction_cell*p->fieldScalar[q]/sum_phi[p->map_sub_to_box[q]]) * vec_v_y[q+store_site_velocity_index[i+other_site_box*store_max_size]];
					V_total_x += (double)(-friction_cell) * vec_v_x[store_site_velocity_index[i+other_site_box*store_max_size]];
					V_total_y += (double)(-friction_cell) * vec_v_y[store_site_velocity_index[i+other_site_box*store_max_size]];
					V_total_x += (double)(friction_cell) * vec_v_x[q+field_start_index[p->index]];
					V_total_y += (double)(friction_cell) * vec_v_y[q+field_start_index[p->index]];
				}
			}
		}
	}
	std::cout<<"velocities: "<<V_total_x<<" "<<V_total_y<<" "<< F_total_x<<" "<<F_total_y <<" "<<V_CoM_x<<" "<<V_CoM_y<<" "<<F_CoM_x<<" "<<F_CoM_y<<" "<<V_All_x<<" "<<V_All_y<<" "<<V_phi_x<<" "<<V_phi_y<<" "<<F_phi_x<<" "<<F_phi_y <<" " << vec_f_x[450]<<" "<< vec_f_y[450] << " "<< V_gammaphi_x<< " "<< V_gammaphi_y  <<std::endl;*/

}

void WetModel::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sum_phi[k]+=p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;

	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));

        /*p->fieldDX[q] = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
	if(box->getWalls(k)<wall_slip){
        	p->fieldDY[q] = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
		p->laplacianPhi[q] = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[7+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + p->fieldScalar[p->neighbors_sub[1+q*9]] - 4.*p->fieldScalar[q];

		//p->Phi00[q] = (p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] - 2 * p->fieldScalar[q]) - (p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] - 2 * p->fieldScalar[q]);
		//p->Phi01[q] = 0.25 * (p->fieldScalar[p->neighbors_sub[8+q*9]] - p->fieldScalar[p->neighbors_sub[6+q*9]] - p->fieldScalar[p->neighbors_sub[2+q*9]] + p->fieldScalar[p->neighbors_sub[0+q*9]]);
	}
	else{
        	p->fieldDY[q] = 0.;
		p->laplacianPhi[q] = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + - 2.*p->fieldScalar[q];

		//p->Phi00[q] = (p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] - 2 * p->fieldScalar[q]);
		//p->Phi01[q] = 0.;
	}*/

        p->fieldDX[q] = BaseInteraction::derivX(p, q, k);
        p->fieldDY[q] = BaseInteraction::derivY(p, q, k);
	p->laplacianPhi[q] = BaseInteraction::Laplacian(p, q, k);

	BaseInteraction::updateFieldProperties(p, q, k);

	//sumQ00[k] += (p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]);
	//sumQ01[k] += (p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]);

	if(size_store_site_velocity_index[k]>=store_max_size){
		for(int m=0; m<size_store_site_velocity_index[k];m++){
			std::cout<<"Too many fields list: "<<store_site_velocity_index[m+k*store_max_size]<<std::endl;
		}	
		throw RCexception("Too many field patches overlap: %d, %d, %d, %d, %d, %d", p->index, k, q, p->sub_corner_bottom_left, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
	}

	store_site_velocity_index[size_store_site_velocity_index[k]+k*store_max_size]=q+field_start_index[p->index];
	store_site_field[size_store_site_velocity_index[k]+k*store_max_size]=p->fieldScalar[q];
	size_store_site_velocity_index[k]++;
}

number WetModel::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
	if(sum_phi[k] - p->fieldScalar[q] > 0){
		p->NQ00 += (p->fieldScalar[q] / p->area) * (sumQ00[k] - p->fieldScalar[q] * p->Q00) / (sum_phi[k] - p->fieldScalar[q]);
		p->NQ01 += (p->fieldScalar[q] / p->area) * (sumQ01[k] - p->fieldScalar[q] * p->Q01) / (sum_phi[k] - p->fieldScalar[q]);
	}

        //number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        //number dy = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );

	//number ytop, ybottom, xleft, xright;
	/*if(p->neighbors_sub[5+q*9]==-1)xright=p->fieldScalar[q];
	else xright=p->fieldScalar[p->neighbors_sub[5+q*9]];
	if(p->neighbors_sub[3+q*9]==-1)xleft=p->fieldScalar[q];
	else xleft=p->fieldScalar[p->neighbors_sub[3+q*9]];
	if(p->neighbors_sub[7+q*9]==-1)ytop=p->fieldScalar[q];
	else ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
	if(p->neighbors_sub[1+q*9]==-1)ybottom=p->fieldScalar[q];
	else ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
	if(p->neighbors_sub[5+q*9]==-1)xright=0;
	else xright=p->fieldScalar[p->neighbors_sub[5+q*9]];
	if(p->neighbors_sub[3+q*9]==-1)xleft=0;
	else xleft=p->fieldScalar[p->neighbors_sub[3+q*9]];
	if(p->neighbors_sub[7+q*9]==-1)ytop=0;
	else ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
	if(p->neighbors_sub[1+q*9]==-1)ybottom=0;
	else ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
	number dx = .5*( xright - xleft );
	number dy = .5*( ytop - ybottom );*/


        //p->fieldDX[q] = dx;
        //p->fieldDY[q] = dy;

	//this part gets the field values in teh respective directions from q;
	//It is hardcoded so take care, the relvant part is that the lattice is square;
	//The neighbors start form the top and rotate couterclockwise.

	//number laplacianPhi = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[7+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + p->fieldScalar[p->neighbors_sub[1+q*9]] - 4.*p->fieldScalar[q];
	//number laplacianPhi = xright + ytop + xleft + ybottom - 4.*p->fieldScalar[q];

	number laplacianSquare;
	if(box->getWalls(k)<wall_slip){
		laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];
	}
	else{
		laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[3+k*9]] - 2.*phi2[k];
	}

	// CH term coupled to chemical (use first)
	number CH = gamma*(8*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-2*p->fieldScalar[q])/lambda - 2*lambda*p->laplacianPhi[q]);
	//number CH = gamma*(8*p->fieldScalar[q]*(p->fieldScalar[q]-1)*(2*p->fieldScalar[q]-1)/lambda - 2*lambda*laplacianPhi);
	//CH term anisotropy
	/*number phixx = (p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] - 2 * p->fieldScalar[q]);
	number phiyy = (p->fieldScalar[p->neighbors_sub[7+q*9]] + p->fieldScalar[p->neighbors_sub[1+q*9]] - 2 * p->fieldScalar[q]);
	number phixy = 0.25 * (p->fieldScalar[p->neighbors_sub[8+q*9]] - p->fieldScalar[p->neighbors_sub[6+q*9]] - p->fieldScalar[p->neighbors_sub[2+q*9]] + p->fieldScalar[p->neighbors_sub[0+q*9]]);
	number norm = sqrt(p->nemQ[0]*p->nemQ[0] + p->nemQ[1]*p->nemQ[1]);
	number v0 = p->nemQ[0] / norm;
	number v1 = p->nemQ[1] / norm;	
	number k00 = R1 * v0 * v0 + R2 * (1 - v0 * v0);
	number k11 = R1 * v1 * v1 + R2 * (1 - v1 * v1);
	number k01 = R1 * v1 * v0 - R2 * v1 * v0;
	number anisotropy = k00 * phixx + 2 * k01 * phixy + k11 * phiyy;
	number CH = gamma*(8*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-2*p->fieldScalar[q])/lambda - lambda * anisotropy);*/
	//if(p->index==0 && q==0)std::cout<<sqrt(p->nemQ[0]*p->nemQ[0] + p->nemQ[1]*p->nemQ[1])<<" "<<sqrt(v0*v0+v1*v1)  <<std::endl;
	
   
	// area conservation term
	number A = - 4*(mu/a0)*(1-p->area/a0)*p->fieldScalar[q];

	// repulsion term
	number Rep = 4*(kappa/lambda)*p->fieldScalar[q]*(phi2[k]-p->fieldScalar[q]*p->fieldScalar[q]);

	// adhesion term
	number lsquare = 2 * p->fieldScalar[q] * p->laplacianPhi[q] + 2 * (p->fieldDX[q] * p->fieldDX[q] + p->fieldDY[q] * p->fieldDY[q]);
	number suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
	number Adh = - 4*lambda*omega*suppress*p->fieldScalar[q];


	//shape term
	/*number x_x_com = p->map_sub_to_box_x[q] - p->CoM[0];
	number y_y_com = p->map_sub_to_box_y[q] - p->CoM[1];
	if(x_x_com >= box->getXsize()/2)x_x_com -= box->getXsize();
	if(x_x_com <= -box->getXsize()/2)x_x_com += box->getXsize();
	if(y_y_com >= box->getYsize()/2)y_y_com -= box->getYsize();
	if(y_y_com <= -box->getYsize()/2)y_y_com += box->getYsize();

	number val1 = (x_x_com * x_x_com) * ( (1 + p->Q00) / R1 + (1 - p->Q00) / R2 );
	number val2 = (2 * x_x_com * y_y_com) * ( p->Q01 / R1 + p->Q01 / R2 );
	number val3 = (y_y_com * y_y_com) * ( (1 - p->Q00) / R1 + (1 + p->Q00) / R2 );

	number Shape = 2 * Kg * (p->fieldScalar[q] - exp( -(val1 + val2 + val3) ) );*/
	//number Shape = 2 * Kg * ( p->fieldScalar[q] - exp( - ((x_x_com * x_x_com / R1 + y_y_com * y_y_com / R2) * (1 + p->Q00) + (x_x_com * x_x_com / R2 + y_y_com * y_y_com / R1) * (1 - p->Q00) + 2 * x_x_com * y_y_com * p->Q01 * R_term) ) );

	/*number D_i = sqrt(p->S00 * p->S00 + p->S01 * p->S01);
	number dDidphi = 2 * p->S00 * p->Phi00[q] + 4 * p->S01 * p->Phi01[q];
	//number Shape = - 2 * Kg * (2 * D_i - D_0) * dDidphi / (D_i * D_i * D_i);
	number Shape;
	if(D_i!=0)Shape = 2 * Kg * (2 * D_i - 3.) * dDidphi / (D_i);
	else Shape = 0.;*/
	//if(p->index==0 && q==0)std::cout<<Kg<<" "<<Shape<<" "<<D_i <<std::endl;



	// delta F / delta phi_i
	number V = CH + A + Rep + Adh;
	//number V = CH + A + Rep + Adh + Shape;
	//if(p->index==0 && q==0)std::cout<<p->freeEnergy[q]<<" "<<p->LsubX<<" "<<p->LsubY<<" "<< p->neighbors_sub[5+q*9]<<" "<< p->neighbors_sub[3+q*9]<<" "<< p->neighbors_sub[7+q*9]<<" "<<p->neighbors_sub[1+q*9]<<" "<<CH<<" "<<A<<" "<<Rep<<std::endl;
	p->freeEnergy[q] += V;
	p->Pressure[q] = Rep - CH - A;
	//if(p->index==0 && q==0)std::cout<<"Free energy: "<<p->freeEnergy[q]<<" "<<CH<<" "<<A<<" "<<Rep<<" "<<Adh<<" "<<dx<<" "<<dy <<std::endl;

	return V;
}


void WetModel::calc_internal_forces(BaseField *p, int q) {

        //int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];

	//if(p->index==0 && q==0)std::cout<<"Forces: "<< p->freeEnergy[q]*p->fieldDX[q] <<" "<< p->freeEnergy[q]*p->fieldDY[q]<<" "<<  0.5 * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] )  << " "<< 0.5 * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] ) <<std::endl;

	//passive (passive force)
	number f_passive_x = (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] );
	number f_passive_y = (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] );
	//number f_passive_x = p->freeEnergy[q]*p->fieldDX[q];
	//number f_passive_y = p->freeEnergy[q]*p->fieldDY[q];
	p->Fpassive_x[q] = f_passive_x * passive_alpha;
	p->Fpassive_y[q] = f_passive_y * passive_alpha;

	//active inter cells (active force)
	number fQ_self_x = - (p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]);
	number fQ_self_y = - (p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]);

	number fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ01[box->neighbors[7+k*9]] - sumQ01[box->neighbors[1+k*9]] ) ) - fQ_self_x;
	number fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( sumQ00[box->neighbors[7+k*9]] - sumQ00[box->neighbors[1+k*9]] ) ) - fQ_self_y;
	//number fQ_inter_x = - ( sumQ00[k] - (p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]));
	//number fQ_inter_y = - ( sumQ01[k] - (p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]));
	//number fQ_inter_x = - ( sumQ00[k] + fQ_self_x);
	//number fQ_inter_y = - ( sumQ01[k] + fQ_self_y);
	//number fQ_inter_x = - (sumQ00[k]*p->fieldDX[q] + sumQ01[k]*p->fieldDY[q]);
	//number fQ_inter_y = - (sumQ01[k]*p->fieldDX[q] - sumQ00[k]*p->fieldDY[q]);
	//number fQ_inter_x = - (p->NQ00 * p->fieldDX[q] + p->NQ01 * p->fieldDY[q]);
	//number fQ_inter_y = - (p->NQ01 * p->fieldDX[q] - p->NQ00 * p->fieldDY[q]);

	p->Factive_x[q] = zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive_y[q] = zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	//p->velocityX[q] = (p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter)/friction;
	//p->velocityY[q] = (p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter)/friction;

	//vec_f_x[q+p->index*p->subSize] = p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
	//vec_f_y[q+p->index*p->subSize] = p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
	/*double v0=0.01;
	if(p->index==0)v0*=1;
	else v0*=-1;*/

	if(box->getWalls(k)<wall_slip){
		//if(p->index==0){
		//F_total_x += p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		//F_total_y += p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
		//F_total_x += p->fieldScalar[q] * grad_free_energy_x[k] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		//F_total_y += p->fieldScalar[q] * grad_free_energy_y[k] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
        	//F_total_x += (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] ) + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
        	//F_total_y += (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] ) + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		//std::cout<< fQ_self_x * zetaQ_self<<" "<<fQ_self_y * zetaQ_self<<" "<< fQ_inter_x * zetaQ_inter  <<" "<< fQ_inter_y * zetaQ_inter  <<std::endl;
		//}
		//F_total_x += p->freeEnergy[q] * p->fieldDX[q];
		//F_total_y += p->freeEnergy[q] * p->fieldDY[q];
		//if(abs(fQ_self_x * zetaQ_self)>F_total_x) F_total_x = fQ_self_x * zetaQ_self;
		//if(abs(fQ_self_y * zetaQ_self)>F_total_y) F_total_y = fQ_self_y * zetaQ_self;
		//F_total_x += fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		//F_total_y += fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
		//F_total_x += f_passive_x + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		//F_total_y += f_passive_y + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;

		//vec_f_x[q+field_start_index[p->index]] = p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		//vec_f_y[q+field_start_index[p->index]] = p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
		vec_f_x[q+field_start_index[p->index]] = f_passive_x * passive_alpha + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		vec_f_y[q+field_start_index[p->index]] = f_passive_y * passive_alpha + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
	}
	else{
		vec_f_x[q+field_start_index[p->index]] = 0.;
		vec_f_y[q+field_start_index[p->index]] = 0.;
	}
}


void WetModel::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	if(store)p->thetaQ_old = p->thetaQ;

	p->thetaQ = p->thetaQ_old - dt * J_Q * atan2(p->S00 * p->Q01 - p->S01 * p->Q00, p->S00 * p->Q00 + p->S01 * p->Q01);
	p->Q00 = cos(2 * p->thetaQ);
	p->Q01 = sin(2 * p->thetaQ);

    	/*number nemQ_mod = sqrt(p->Q00 * p->Q00 + p->Q01 * p->Q01);
    	p->nemQ[0] = sqrt((1 + p->Q00/nemQ_mod)/2);
	number sgn;
	if(p->Q01>0)sgn=1;
	else if(p->Q01<0) sgn=-1;
	else sgn=0;
    	p->nemQ[1] = sgn*sqrt((1 - p->Q00/nemQ_mod)/2);*/


    	/*number S_mod = sqrt(p->S00 * p->S00 + p->S01 * p->S01);
	if(store)p->nemQ_old = {p->nemQ[0] , p->nemQ[1]};
	if(S_mod>0.00000001){
		number t = 0.5 * atan2(p->S01, p->S00);
		std::vector<number> d = {cos(t) , sin(t)};
		number sgn = (d[0] * p->nemQ[0] + d[1] * p->nemQ[1] > 0.0)? 1.0:-1.0;
		p->nemQ[0] = p->nemQ_old[0] + dt * J_Q * (sgn * d[0] - p->nemQ[0]);
		p->nemQ[1] = p->nemQ_old[1] + dt * J_Q * (sgn * d[1] - p->nemQ[1]);

		p->Q00 = 0.5 * (p->nemQ[0] * p->nemQ[0] - p->nemQ[1] * p->nemQ[1]);
		p->Q01 = p->nemQ[0] * p->nemQ[1];
	}*/

	/*p->Q00 += dt * J_Q * (p->S00 - p->Q00);
	p->Q01 += dt * J_Q * (p->S01 - p->Q01);
    	number nemQ_mod = sqrt(p->Q00 * p->Q00 + p->Q01 * p->Q01);
	if(nemQ_mod>0.000000001){
	    	number nx = sqrt((1 + p->Q00/nemQ_mod)/2);
		number sgn;
		if(p->Q01>0)sgn=1;
		else if(p->Q01<0) sgn=-1;
		else sgn=0;
    		number ny = sgn*sqrt((1 - p->Q00/nemQ_mod)/2);
		p->nemQ[0]=nemQ_mod * nx;
		p->nemQ[1]=nemQ_mod * ny;
	}
	else{
		p->nemQ[0]=0.;
		p->nemQ[1]=0.;
	}*/

	if(anchoring) update_anchoring(p);
}


void WetModel::update_anchoring(BaseField*p){

	number walls_length = 8;
	number dist1 = ((double)p->CoM_old[1] - walls_length);
	number dist2 = ((double)box->getYsize() - walls_length) - (double)p->CoM_old[1];
	number theta;
	if(dist1<dist2) theta = (PI/2) - (PI/2) * dist1 / (((double)box->getYsize() - 2 * walls_length)/2);
	else theta = (PI/2) + ((PI/2) * dist2 / (((double)box->getYsize() - 2 * walls_length)/2));

	//std::cout<<"values: "<<p->index<<" "<<p->CoM[1]<<" "<<theta<<std::endl;

	/*p->Q00=0.5*cos(2*theta);
	p->Q01=0.5*sin(2*theta);
    	number nemQ_mod = sqrt(p->Q00 * p->Q00 + p->Q01 * p->Q01);
    	number nx = sqrt((1 + p->Q00/nemQ_mod)/2);
	number sgn;
	if(p->Q01>0)sgn=1;
	else if(p->Q01<0) sgn=-1;
	else sgn=0;
    	number ny = sgn*sqrt((1 - p->Q00/nemQ_mod)/2);
	p->nemQ[0]=nx;
	p->nemQ[1]=ny;*/

	p->nemQ[0]=cos(theta);
	p->nemQ[1]=sin(theta);
	p->Q00 = 0.5 * (p->nemQ[0] * p->nemQ[0] - p->nemQ[1] * p->nemQ[1]);
	p->Q01 = p->nemQ[0] * p->nemQ[1];


	/*number delta = -PI/12;
	number theta = 0;
	number S = 0.5;

	if(p->CoM[1] < 2.5 * R){
		theta = (PI/2 + delta);
		p->Q00 = S * cos(2*theta);
		p->Q01 = S * sin(2*theta);
	}
	else if(p->CoM[1] > box->getYsize() - 2.5 * R){
		theta = -(PI/2 + delta);
		p->Q00 = S * cos(2*theta);
		p->Q01 = S * sin(2*theta);
	}*/
}
