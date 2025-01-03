#include "WetModel.h"
#include "../Utilities/Utils.h"

WetModel::WetModel() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				omega(0.004),
				mu(3.),
				kappa(0.1),
				friction(1.),
				zetaQ_self(0),
				zetaQ_inter(0),
				J_Q(1),
				friction_cell(1.),
				tolerance(0.0001),
				wall_slip(0.5){
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
	getInputNumber(&inp, "J_Q", &J_Q, 0);
	getInputBool(&inp, "anchoring", &anchoring, 0);
	getInputNumber(&inp, "friction_cell", &friction_cell_active, 0);
	getInputNumber(&inp, "friction", &friction_active, 0);
	getInputNumber(&inp, "CGtolerance", &tolerance, 0);
	getInputNumber(&inp, "wall_slip", &wall_slip, 0);
}

void WetModel::init() {
        a0=PI*R*R;
	store_max_size=20;
	//solverCG.setTolerance(0.000000001);
	if(tolerance>0)
		solverCG.setTolerance(tolerance);
        set_omp_tasks(omp_thread_num);
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
	zetaQ_self=zetaQ_self_active;
	zetaQ_inter=zetaQ_inter_active;
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
}

void WetModel::initFieldProperties(BaseField *p) {

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


void WetModel::updateFieldProperties(BaseField *p, int q, int k) {
	BaseInteraction::updateFieldProperties(p, q, k);
	number dx = p->fieldDX[q]; 
	number dy = p->fieldDY[q]; 
	p->S00 += -0.5*(dx*dx-dy*dy);
	p->S01 += -dx*dy;
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
	size_rows_old = size_rows;
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
                p->Factive = std::vector<number> {0., 0.};
                p->Fpassive = std::vector<number> {0., 0.};
                for(int q=0; q<p->subSize;q++){
			calc_internal_forces(p, q);

			//populate sparse matrix
			if(box->getWalls(p->map_sub_to_box[q])<wall_slip)
				//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction+4*friction_cell)));
				tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction+8*friction_cell)));
				//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction)));
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
				other_site_box = box->neighbors[j+p->map_sub_to_box[q]*9];
				//if(sum_phi[q]==0)continue;

				//if(box->getWalls(p->map_sub_to_box[q])<wall_slip)
					//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], other_site_patch+field_start_index[p->index] , (double)(friction_cell)));

				for(int i=0; i<size_store_site_velocity_index[other_site_box]; i++){
					//number weight = p->fieldScalar[q] * store_site_field[i+other_site_box*store_max_size]/(sum_phi[p->map_sub_to_box[q]] * sum_phi[other_site_box]);
					//tri_t_x.push_back(Eigen::Triplet<double> (store_site_velocity_index[i+other_site_box*store_max_size], q+field_start_index[p->index], (double)(-friction_cell*p->fieldScalar[other_site_patch]/sum_phi[other_site_box])));
					
					if(box->getWalls(p->map_sub_to_box[q])<wall_slip){
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], q+field_start_index[p->index], (double)(friction_cell)*weight));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell)*weight));
						tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell*store_site_field[i+other_site_box*store_max_size]/sum_phi[other_site_box])));
						//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(-friction_cell*p->fieldScalar[q]/sum_phi[p->map_sub_to_box[q]])));
					}

					//tri_t_x.push_back(Eigen::Triplet<double> (q+field_start_index[p->index], store_site_velocity_index[i+other_site_box*store_max_size], (double)(friction_cell*p->fieldScalar[other_site_patch]/sum_phi[other_site_box])));
				}
			}
		}
        }

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

	//std::cout<<Eigen::nbThreads( )<<std::endl;
	//std::cout<<"start solver: "<<size_rows/2 <<std::endl;
	mat_m_x.setFromTriplets(tri_t_x.begin(), tri_t_x.end());
	solverCG.compute(mat_m_x);
	//if(size_rows != size_rows_old){
	vec_v_x = solverCG.solve(vec_f_x);
	vec_v_y = solverCG.solve(vec_f_y);
	//}
	/*else{
		vec_v_x = solverCG.solveWithGuess(vec_f_x, vec_v_x);
		vec_v_y = solverCG.solveWithGuess(vec_f_y, vec_v_y);
	}*/
	//std::cout << "#iterations:     " << solverCG.iterations() << std::endl;
	//std::cout << "estimated error: " << solverCG.error()      << std::endl;	
	//std::cout<<"end solver"<<std::endl;
}

void WetModel::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sum_phi[k]+=p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;

	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));

	if(size_store_site_velocity_index[k]>=store_max_size){
		for(int m=0; m<size_store_site_velocity_index[k];m++){
			std::cout<<"Too many fields list: "<<store_site_velocity_index[m+k*store_max_size]<<std::endl;
		}	
		throw RCexception("Too many field patches overlap: %d, %d, %d, %d, %d, %d", p->index, k, q, p->sub_corner_bottom_left, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
	}
	//store_site_velocity_index[size_store_site_velocity_index[k]+k*store_max_size]=q+p->index*p->subSize;
	store_site_velocity_index[size_store_site_velocity_index[k]+k*store_max_size]=q+field_start_index[p->index];
	store_site_field[size_store_site_velocity_index[k]+k*store_max_size]=p->fieldScalar[q];
	size_store_site_velocity_index[k]++;
}

number WetModel::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
        p->fieldDX[q] = dx;
        p->fieldDY[q] = dy;

	//this part gets the field values in teh respective directions from q;
	//It is hardcoded so take care, the relvant part is that the lattice is square;
	//The neighbors start form the top and rotate couterclockwise.
	number laplacianPhi = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[7+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + p->fieldScalar[p->neighbors_sub[1+q*9]] - 4.*p->fieldScalar[q];
	number laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];

	// CH term coupled to chemical
	number CH =+ gamma*(8*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-2*p->fieldScalar[q])/lambda - 2*lambda*laplacianPhi);
   
	// area conservation term
	number A = - 4*mu/a0*(1-p->area/a0)*p->fieldScalar[q];

	// repulsion term
	number Rep = + 4*kappa/lambda*p->fieldScalar[q]*(phi2[k]-p->fieldScalar[q]*p->fieldScalar[q]);

	// adhesion term
	number lsquare = 2 * p->fieldScalar[q] * laplacianPhi + 2 * (dx * dx + dy * dy);
	number suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
	number Adh = - 4*lambda*omega*suppress*p->fieldScalar[q];

	// delta F / delta phi_i
	number V = CH + A + Rep + Adh;
	p->freeEnergy[q] += V;

	return V;
}


void WetModel::calc_internal_forces(BaseField *p, int q) {

        //int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];

	//passive (passive force)
	p->Fpassive[0] += p->freeEnergy[q]*p->fieldDX[q];
	p->Fpassive[1] += p->freeEnergy[q]*p->fieldDY[q];

	//active inter cells (active force)
	number fQ_self_x = -(p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]);
	number fQ_self_y = -(p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]);

	number fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ01[box->neighbors[7+k*9]] - sumQ01[box->neighbors[1+k*9]] ) ) - fQ_self_x;
	number fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( sumQ00[box->neighbors[7+k*9]] - sumQ00[box->neighbors[1+k*9]] ) ) - fQ_self_y;

	p->Factive[0] += zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive[1] += zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	//p->velocityX[q] = (p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter)/friction;
	//p->velocityY[q] = (p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter)/friction;

	//vec_f_x[q+p->index*p->subSize] = p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
	//vec_f_y[q+p->index*p->subSize] = p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
	/*double v0=0.01;
	if(p->index==0)v0*=1;
	else v0*=-1;*/

	if(box->getWalls(k)<wall_slip){
		vec_f_x[q+field_start_index[p->index]] = p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
		vec_f_y[q+field_start_index[p->index]] = p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
	}
	else{
		vec_f_x[q+field_start_index[p->index]] = 0.;
		vec_f_y[q+field_start_index[p->index]] = 0.;
	}
}


void WetModel::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	if(store)p->nemQ_old = {p->nemQ[0] , p->nemQ[1]};
	
	number t = 0.5 * atan2(p->S01, p->S00);
	std::vector<number> d = {cos(t) , sin(t)};
	number sgn = (d[0] * p->nemQ[0] + d[1] * p->nemQ[1] > 0.0)? 1.0:-1.0;
	p->nemQ[0] = p->nemQ_old[0] + dt * J_Q * (sgn * d[0] - p->nemQ[0]);
	p->nemQ[1] = p->nemQ_old[1] + dt * J_Q * (sgn * d[1] - p->nemQ[1]);

	p->Q00 = 0.5 * (p->nemQ[0] * p->nemQ[0] - p->nemQ[1] * p->nemQ[1]);
	p->Q01 = p->nemQ[0] * p->nemQ[1];

	if(anchoring) update_anchoring(p);
}


void WetModel::update_anchoring(BaseField*p){

	number delta = -PI/12;
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
	}
}
