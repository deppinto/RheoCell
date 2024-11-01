#include "WetModel.h"
#include "../Utilities/Utils.h"

WetModel::WetModel() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				omega(0.004),
				mu(3.),
				R(8),
				kappa(0.1),
				friction(2.),
				zetaQ_self(0),
				zetaQ_inter(0),
				J_Q(1),
				friction_cell(3.0){
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
	getInputNumber(&inp, "friction", &friction, 0);
	getInputNumber(&inp, "zetaQ_self", &zetaQ_self_active, 0);
	getInputNumber(&inp, "zetaQ_inter", &zetaQ_inter_active, 0);
	getInputNumber(&inp, "J_Q", &J_Q, 0);
	getInputBool(&inp, "anchoring", &anchoring, 0);
	getInputNumber(&inp, "friction_cell", &friction_cell, 0);
}

void WetModel::init() {
        a0=PI*R*R;
	store_max_size=10;
	solverCG.setTolerance(0.000000001);
}

void WetModel::read_topology(std::vector<BaseField*> &fields) {
        int N = fields.size();

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
}

void WetModel::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
	phi2.resize(Lx*Ly);
	sumQ00.resize(Lx*Ly);
	sumQ01.resize(Lx*Ly);
	sum_phi.resize(Lx*Ly);
	size_store_site_velocity_index.resize(Lx*Ly);
	store_site_velocity_index.resize(Lx*Ly*store_max_size);
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
	        number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
	        p->fieldDX[q] = dx;
	        p->fieldDY[q] = dy;

		p->S00 += -0.5*(dx*dx-dy*dy);
		p->S01 += -dx*dy;
	}
	size_rows += p->subSize;
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

	if(size_rows/CONFIG_INFO->fields()[0]->subSize!=2*CONFIG_INFO->N())throw RCexception("Not all fields have the same patch size for the wet model. Aborting");
	vec_v_x.resize(size_rows);
	vec_f_x.resize(size_rows);
	vec_v_y.resize(size_rows);
	vec_f_y.resize(size_rows);
	mat_m_x.resize(size_rows, size_rows);
}

void WetModel::begin_energy_computation(std::vector<BaseField *> &fields) {

	std::fill(size_store_site_velocity_index.begin(), size_store_site_velocity_index.end(), 0);
	for(auto p : fields) {
		for(int q=0; q<p->subSize;q++)
			computeGlobalSums(p, q, false);
	}

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
			tri_t_x.push_back(Eigen::Triplet<double> (q+p->index*p->subSize, q+p->index*p->subSize, (double)(friction+4*friction_cell)));

			for(auto j : neigh_values){
				other_site_patch = p->neighbors_sub[j+q*9];
				other_site_box = p->map_sub_to_box[other_site_patch];
				if(sum_phi[other_site_box]==0)continue;
				for(int i=0; i<size_store_site_velocity_index[other_site_box]; i++){
					index = store_site_velocity_index[i+other_site_box*store_max_size]/p->subSize;
					sub_q = (store_site_velocity_index[i+other_site_box*store_max_size]-(index*p->subSize));
					tri_t_x.push_back(Eigen::Triplet<double> (q+p->index*p->subSize, sub_q+index*p->subSize, (double)(-friction_cell*p->fieldScalar[other_site_patch]/sum_phi[other_site_box])));
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

	//std::cout<<"start solver: "<<size_rows/2 <<std::endl;
	mat_m_x.setFromTriplets(tri_t_x.begin(), tri_t_x.end());
	//mat_m_x.makeCompressed();
	//std::cout<<"copy done"<<std::endl;
	//solver.analyzePattern(mat_m_x);
	//solver.factorize(mat_m_x);
	solverCG.compute(mat_m_x);
	vec_v_x = solverCG.solveWithGuess(vec_f_x, vec_v_x);
	vec_v_y = solverCG.solveWithGuess(vec_f_y, vec_v_y);
	//std::cout << "#iterations:     " << solverCG.iterations() << std::endl;
	//std::cout << "estimated error: " << solverCG.error()      << std::endl;	
	//std::cout<<"end solver"<<std::endl;


	//std::cout<<"testing xy: "<<vec_v[312*2]<<" "<<vec_v[312*2+1]<<std::endl;
	//std::cout<<"testing x: "<<vec_v_x[312]<<" "<<vec_v_y[312]<<std::endl;
	//std::cout<<"testing old: "<< fields[0]->velocityX[312] <<" "<< fields[0]->velocityY[312]<<std::endl;
	//exit (911);

	//velX = p->Fpassive[0] + p->Factive[0];
	//velY = p->Fpassive[1] + p->Factive[1];
	//K += .5 * (velX * velX + velY * velY);
}

void WetModel::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sum_phi[k]+=p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;

	store_site_velocity_index[size_store_site_velocity_index[k]+k*store_max_size]=q+p->index*p->subSize;
	size_store_site_velocity_index[k]++;
	if(size_store_site_velocity_index[k]>=store_max_size)throw RCexception("Too many field patches overlap!");

	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
}

number WetModel::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
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
	number lsquare = 2 * p->fieldScalar[q] * laplacianPhi + 2 * (dx *dx + dy * dy);
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

	number fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ01[box->neighbors[1+k*9]] - sumQ01[box->neighbors[7+k*9]] ) ) - fQ_self_x;
	number fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( sumQ00[box->neighbors[1+k*9]] - sumQ00[box->neighbors[7+k*9]] ) ) - fQ_self_y;

	p->Factive[0] += zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive[1] += zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	//p->velocityX[q] = (p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter)/friction;
	//p->velocityY[q] = (p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter)/friction;

	//vec_f[0+(q+p->index*p->subSize)*2] = p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
	//vec_f[1+(q+p->index*p->subSize)*2] = p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;

	vec_f_x[q+p->index*p->subSize] = p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter;
	vec_f_y[q+p->index*p->subSize] = p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter;
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
