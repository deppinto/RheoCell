#include "WetPolarModel.h"
#include "../Utilities/Utils.h"

WetPolarModel::WetPolarModel() :
				BaseInteraction(),
				gamma(0.06),
				lambda(2.5),
				omega(0.004),
				mu(20.),
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

WetPolarModel::~WetPolarModel() {

}


void WetPolarModel::get_settings(input_file &inp) {
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
	getInputNumber(&inp, "passive_alpha", &passive_alpha, 0);
}

void WetPolarModel::init() {
        a0=PI*R*R;
	store_max_size=20;
	if(tolerance>0)
		solverCG.setTolerance(tolerance);
        set_omp_tasks(omp_thread_num);
	restart_solver=0;
	std::cout<<"TESTING: Running with tolerance (CG): "<<tolerance<<std::endl;
}

void WetPolarModel::read_topology(std::vector<BaseField*> &fields) {
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
		fields[i]->clock = int(2 * drand48() * J_Q_active);
        }
}

void WetPolarModel::allocate_fields(std::vector<BaseField *> &fields) {
        for(int i = 0; i < (int) fields.size(); i++) {
                fields[i] = new MultiPhaseField();
        }
}

void WetPolarModel::apply_changes_after_equilibration(){

	BaseInteraction::apply_changes_after_equilibration();
	zetaQ_self=zetaQ_self_active;
	zetaQ_inter=zetaQ_inter_active;
	J_Q=J_Q_active;
	friction=friction_active;
	friction_cell=friction_cell_active;
}

void WetPolarModel::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
	if(box->lees_edwards)throw RCexception("Interaction is not compatible with LEBc. Aborting");
	phi2.resize(Lx*Ly);
	sum_phi.resize(Lx*Ly);
	size_store_site_velocity_index.resize(Lx*Ly);
	store_site_velocity_index.resize(Lx*Ly*store_max_size);
	store_site_field.resize(Lx*Ly*store_max_size);
	for(int i =0; i<Lx*Ly; i++){resetSums(i);size_store_site_velocity_index[i]=0;}
}

void WetPolarModel::resetSums(int k) {
	phi2[k]=0;
        sum_phi[k]=0;
}

void WetPolarModel::initFieldProperties(BaseField *p) {

	for(int q=0; q<p->subSize;q++) {
		int k = p->GetSubIndex(q, box);
		BaseInteraction::updateFieldProperties(p, q, k);

	        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
	        number dy = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
	        p->fieldDX[q] = dx;
	        p->fieldDY[q] = dy;

		p->S00 += -0.5*(dx*dx-dy*dy);
		p->S01 += -dx*dy;
	}
}


void WetPolarModel::updateFieldProperties(BaseField *p, int q, int k) {

	p->S00 += -0.5*(p->fieldDX[q]*p->fieldDX[q]-p->fieldDY[q]*p->fieldDY[q]);
	p->S01 += -p->fieldDX[q]*p->fieldDY[q];
}


void WetPolarModel::check_input_sanity(std::vector<BaseField *> &fields) {

}


void WetPolarModel::begin_energy_computation() {
		
        for(int i = 0; i < CONFIG_INFO->N(); i++) {
                initFieldProperties(CONFIG_INFO->fields()[i]);
        }
}

void WetPolarModel::begin_energy_computation(std::vector<BaseField *> &fields) {

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
		std::cout<<p->index<<" "<<p->perimeter<<std::endl;
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


        U = (number) 0;
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++)
			U += f_interaction(p, q);

        }


        K = (number) 0;
	std::vector<Eigen::Triplet<double>> tri_t_x;
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
			if(box->getWalls(p->map_sub_to_box[q])<wall_slip)
				tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction)));
			else
				tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], 1.0));

			for(auto j : neigh_values){
				other_site_box = box->neighbors[j+p->map_sub_to_box[q]*9];

				for(int i=0; i<size_store_site_velocity_index[other_site_box]; i++){

					if(box->getWalls(p->map_sub_to_box[q])<wall_slip){
						tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction_cell * weight_values[j])));
						tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell * weight_values[j])));
					}
				}
			}
		}
        }

	if(friction_cell!=0){
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
}

void WetPolarModel::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sum_phi[k]+=p->fieldScalar[q];
	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));

        p->fieldDX[q] = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
	if(box->getWalls(k)<wall_slip){
        	p->fieldDY[q] = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
		p->laplacianPhi[q] = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[7+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + p->fieldScalar[p->neighbors_sub[1+q*9]] - 4.*p->fieldScalar[q];
	}
	else{
        	p->fieldDY[q] = 0.;
		p->laplacianPhi[q] = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + - 2.*p->fieldScalar[q];
	}
	p->perimeter += p->fieldDX[q] * p->fieldDX[q] + p->fieldDY[q] * p->fieldDY[q];

	BaseInteraction::updateFieldProperties(p, q, k);

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

number WetPolarModel::f_interaction(BaseField *p, int q) {

	int k = p->map_sub_to_box[q];

	number laplacianSquare;
	if(box->getWalls(k)<wall_slip){
		laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];
	}
	else{
		laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[3+k*9]] - 2.*phi2[k];
	}

	// CH term coupled to chemical (use first)
	number CH = gamma*(8*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-2*p->fieldScalar[q])/lambda - 2*lambda*p->laplacianPhi[q]);
   
	// area conservation term
	number A = - 4*(mu/a0)*(1-p->area/a0)*p->fieldScalar[q];

	// repulsion term
	number Rep = 4*(kappa/lambda)*p->fieldScalar[q]*(phi2[k]-p->fieldScalar[q]*p->fieldScalar[q]);

	// adhesion term
	number lsquare = 2 * p->fieldScalar[q] * p->laplacianPhi[q] + 2 * (p->fieldDX[q] * p->fieldDX[q] + p->fieldDY[q] * p->fieldDY[q]);
	number suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
	number Adh = - 4*lambda*omega*suppress*p->fieldScalar[q];

	//Perimeter term
	number p0 = 12; //2 * PI * 8;
	number beta = mu * 10.;
	number P = 4*(beta/p0) * (1 - p->perimeter/p0) * p->laplacianPhi[q];
	//number P = 0;
	

	// delta F / delta phi_i
	number V = CH + A + Rep + Adh + P;
	p->freeEnergy[q] += V;
	p->Pressure[q] = Rep - CH - A;
	return V;
}


void WetPolarModel::calc_internal_forces(BaseField *p, int q) {

	int k = p->map_sub_to_box[q];

	//passive (passive force)
	number f_passive_x = (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] );
	number f_passive_y = (-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] );
	p->Fpassive_x[q] = f_passive_x * passive_alpha;
	p->Fpassive_y[q] = f_passive_y * passive_alpha;

	//active inter cells (active force)
	double Theta =  (-p->fieldDX[q] * p->Q00) < 0 ? 0.0 : 1.0;
	number fQ_self_x = lambda * (-p->fieldDX[q] * p->Q00) * Theta * p->Q00;
	Theta =  (-p->fieldDY[q] * p->Q01) < 0 ? 0.0 : 1.0;
	number fQ_self_y = lambda * (-p->fieldDY[q] * p->Q01) * Theta * p->Q01;

	p->Factive_x[q] = zetaQ_self * fQ_self_x;
	p->Factive_y[q] = zetaQ_self * fQ_self_y;


	if(box->getWalls(k)<wall_slip){
		vec_f_x[q+field_start_index[p->index]] = f_passive_x * passive_alpha + fQ_self_x * zetaQ_self;
		vec_f_y[q+field_start_index[p->index]] = f_passive_y * passive_alpha + fQ_self_y * zetaQ_self;
	}
	else{
		vec_f_x[q+field_start_index[p->index]] = 0.;
		vec_f_y[q+field_start_index[p->index]] = 0.;
	}
}


void WetPolarModel::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	//if(store)p->thetaQ_old = p->thetaQ;
	//p->thetaQ = p->thetaQ_old + sqrt(dt) * J_Q * Utils::gaussian();
	//p->Q00 = cos(p->thetaQ);
	//p->Q01 = sin(p->thetaQ);

	p->clock += 1;
	//std::cout<<p->clock<<std::endl;
	if(p->clock == J_Q_active * 2)
	{
		//std::cout<<p->index<<" "<<p->clock<<std::endl;
		p->thetaQ = PI * (1-2*drand48());
		p->Q00 = cos(p->thetaQ);
		p->Q01 = sin(p->thetaQ);
		p->clock = 0;
	}


	if(anchoring) update_anchoring(p);
}


void WetPolarModel::update_anchoring(BaseField*p){

	number walls_length = 8;
	number dist1 = ((double)p->CoM_old[1] - walls_length);
	number dist2 = ((double)box->getYsize() - walls_length) - (double)p->CoM_old[1];
	number theta;
	if(dist1<dist2) theta = (PI/2) - (PI/2) * dist1 / (((double)box->getYsize() - 2 * walls_length)/2);
	else theta = (PI/2) + ((PI/2) * dist2 / (((double)box->getYsize() - 2 * walls_length)/2));

	p->nemQ[0]=cos(theta);
	p->nemQ[1]=sin(theta);
	p->Q00 = 0.5 * (p->nemQ[0] * p->nemQ[0] - p->nemQ[1] * p->nemQ[1]);
	p->Q01 = p->nemQ[0] * p->nemQ[1];
}
