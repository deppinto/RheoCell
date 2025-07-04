#include "DifferentialAdhesion.h"
#include "../Utilities/Utils.h"

DifferentialAdhesion::DifferentialAdhesion() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				mu(3.),
				kappa(0.1),
				omega(0.),
				zetaQ_self(0),
				zetaQ_inter(0),
				strain_rate(0.) {
	a0=PI*R*R;
}

DifferentialAdhesion::~DifferentialAdhesion() {

}


void DifferentialAdhesion::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "R", &R, 0);
	getInputNumber(&inp, "gamma", &gamma, 0);
	getInputNumber(&inp, "mu", &mu, 0);
	getInputNumber(&inp, "kappa", &kappa, 0);
	getInputNumber(&inp, "omega", &omega, 0);
	getInputNumber(&inp, "zetaQ_self", &zetaQ_self_active, 0);
	getInputNumber(&inp, "zetaQ_inter", &zetaQ_inter_active, 0);
	getInputNumber(&inp, "lees_edwards_shear_rate", &strain_rate_active, 1);
}


void DifferentialAdhesion::read_topology(std::vector<BaseField *> &fields) {
        int N = fields.size();
	//field_start_index.resize(N);

        std::ifstream topology(topology_filename, std::ios::in);
        if(!topology.good()) {
                throw RCexception("Can't read topology file '%s'. Aborting", topology_filename);
        }

        allocate_fields(fields);
	size_rows=0;//EXTRA CARE!!! THIS SHOULDNT BE HARD CODED!!!
        for(int i = 0; i < N; i++) {
                fields[i]->index = i;
		fields[i]->get_interaction_values(R);
		size_rows += 30 * 30;
        }
	//size_rows_old = size_rows;
	//vec_omega.resize(size_rows);
	//std::fill(vec_omega.begin(), vec_omega.end(), -1.);
        //phi_omega.resize(size_rows);
	//std::fill(phi_omega.begin(), phi_omega.end(), -1.);
	//construct_omega(fields);
}


void DifferentialAdhesion::set_box(BaseBox *boxArg) {
                box = boxArg;
		int Lx=box->getXsize();
	        int Ly=box->getYsize();
		//size_store_site_omega_index.resize(Lx*Ly);
		//store_site_omega_index.resize(Lx*Ly*store_max_size);
		//store_site_omega_sub.resize(Lx*Ly*store_max_size);
		//store_site_field.resize(Lx*Ly*store_max_size);
		if(box->lees_edwards)throw RCexception("Interaction is not compatible with LEBc. Aborting");
        	phi2.resize(Lx*Ly);
		sumQ00.resize(Lx*Ly);
		sumQ01.resize(Lx*Ly);
		for(int i =0; i<Lx*Ly; i++){resetSums(i);}//size_store_site_omega_index[i]=0;}
}

void DifferentialAdhesion::init() {
	a0=PI*R*R;
	//store_max_size=20;
}

void DifferentialAdhesion::allocate_fields(std::vector<BaseField *> &fields) {
	for(int i = 0; i < (int) fields.size(); i++) {
		fields[i] = new MultiPhaseField();
	}
}

void DifferentialAdhesion::check_input_sanity(std::vector<BaseField *> &fields) {

}

void DifferentialAdhesion::apply_changes_after_equilibration(){
	strain_rate = strain_rate_active;
	zetaQ_self=zetaQ_self_active;
	zetaQ_inter=zetaQ_inter_active;
	//size_rows_old = 0;
}

void DifferentialAdhesion::resetSums(int k) {
	phi2[k]=0;
        sumQ00[k]=0;
        sumQ01[k]=0;
	//size_store_site_omega_index[k]=0;
}

void DifferentialAdhesion::updateFieldProperties(BaseField *p, int q, int k) {
	BaseInteraction::updateFieldProperties(p, q, k);
	p->S00 += -0.5*(p->fieldDX[q]*p->fieldDX[q]-p->fieldDY[q]*p->fieldDY[q]);
	p->S01 += -p->fieldDX[q]*p->fieldDY[q];
}

void DifferentialAdhesion::begin_energy_computation() {
		
        for(int i = 0; i < CONFIG_INFO->N(); i++) {
                initFieldProperties(CONFIG_INFO->fields()[i]);
        }
}

void DifferentialAdhesion::initFieldProperties(BaseField *p) {

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

void DifferentialAdhesion::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->S00;
        sumQ01[k]+=p->fieldScalar[q]*p->S01;
	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));


        p->fieldDX[q] = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        p->fieldDY[q] = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
	BaseInteraction::updateFieldProperties(p, q, k);

	/*if(size_rows != size_rows_old){
		if(size_store_site_omega_index[k]>=store_max_size){
			for(int m=0; m<size_store_site_omega_index[k];m++){
				std::cout<<"Too many fields list: "<<store_site_omega_index[m+k*store_max_size]<<std::endl;
			}	
			throw RCexception("Too many field patches overlap: %d, %d, %d, %d, %d, %d", p->index, k, q, p->sub_corner_bottom_left, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
		}
		store_site_omega_index[size_store_site_omega_index[k]+k*store_max_size] = p->index;
		store_site_omega_sub[size_store_site_omega_index[k]+k*store_max_size] = q;
		store_site_field[size_store_site_omega_index[k]+k*store_max_size] = box->sqr_min_image_distance(std::vector<number>{(number)(k - int(k/box->getXsize()) * box->getXsize()), (number)(int(k/box->getXsize()))}, std::vector<number> {p->CoM[0], p->CoM[1]});
		size_store_site_omega_index[k]++;
	}*/
}

void DifferentialAdhesion::construct_omega(std::vector<BaseField *> &fields) {

	for(auto p : fields) {
		for(int q=0; q<p->subSize;q++){
			int k = p->map_sub_to_box[q];
			//number value_max = -1.;
			number value_min = 10000000;
			bool found_neigh = false;
			for(int z=0; z<size_store_site_omega_index[k]; z++){
				/*if(store_site_field[z + k * store_max_size] > value_max && store_site_omega_index[z + k * store_max_size]!=p->index){
					vec_omega[q + field_start_index[p->index]] = store_site_omega_index[z + k * store_max_size];
					phi_omega[q + field_start_index[p->index]] = store_site_omega_sub[z + k * store_max_size];
				 	value_max = store_site_field[z + k * store_max_size];
					found_neigh = true;
				}*/
				if(store_site_omega_index[z + k * store_max_size] != p->index){
					if(store_site_field[z + k * store_max_size] < value_min){
						vec_omega[q + field_start_index[p->index]] = store_site_omega_index[z + k * store_max_size];
						phi_omega[q + field_start_index[p->index]] = store_site_omega_sub[z + k * store_max_size];
						value_min = store_site_field[z + k * store_max_size];
						found_neigh = true;
					}
				}
			}
			if(found_neigh == false){
				vec_omega[q + field_start_index[p->index]] = 0.;
				phi_omega[q + field_start_index[p->index]] = -1.;
			}
			//else if(p->index==55){std::cout<< q << " "<<vec_omega[q + field_start_index[p->index]] <<std::endl;}
		}
	}
}

void DifferentialAdhesion::begin_energy_computation(std::vector<BaseField *> &fields) {


	for(auto p : fields) {
		for(int q=0; q<p->subSize;q++){
			computeGlobalSums(p, q, false);
		}
	}


	/*if(size_rows != size_rows_old){
		vec_omega.resize(size_rows);
        	phi_omega.resize(size_rows);
		construct_omega(fields);
	}*/
	//std::fill(vec_omega.begin(), vec_omega.end(), 0.);
	//std::fill(phi_omega.begin(), phi2_omega.end(), 0.);
	//size_rows_old = size_rows;
	//size_rows=0;


        U = (number) 0;
        for(auto p : fields) {
		number velocity_value = 0.;
		if(p->index % 10 == 0) velocity_value = -1.;
		if((p->index + 1) % 10 == 0) velocity_value = 1.;
                for(int q=0; q<p->subSize;q++) {
			p->velocityX[q] = 0.;
			p->velocityY[q] = 0.;
			BaseField *pp = fields[0]; //fields[vec_omega[q + field_start_index[p->index]]];
			int qq = 0; //phi_omega[q + field_start_index[p->index]];
                        //U += f_interaction(p, q);
			if(velocity_value < -0.5 || velocity_value > 0.5){
				//int x = q - int(q/p->LsubX) * p->LsubX;			
				if(velocity_value < -0.5)p->velocityX[q] = -1. * strain_rate;
				else if(velocity_value > 0.5)p->velocityX[q] = 1. * strain_rate;
				else U += f_interaction(p, q, pp, qq);
			}
			else U += f_interaction(p, q, pp, qq);
                }
        }

        K = 0.;
	number velX = 0.;
	number velY = 0.;
        for(auto p : fields) {
		number velocity_value = 0.;
		if(p->index % 10 == 0) velocity_value = -1.;
		if((p->index + 1) % 10 == 0) velocity_value = 1.;
                for(int q=0; q<p->subSize;q++){
			if(velocity_value < -0.5 || velocity_value > 0.5){
				if(velocity_value < -0.5)p->velocityX[q] = -1. * strain_rate;
				else if(velocity_value > 0.5)p->velocityX[q] = 1. * strain_rate;
				else calc_internal_forces(p, q);
			}
			else calc_internal_forces(p, q);

			//if(p->index % 10 != 0 && (p->index + 1) % 10 != 0)calc_internal_forces(p, q);
                	velX += p->Fpassive_x[q] + p->Factive_x[q];
                	velY += p->Fpassive_y[q] + p->Factive_y[q];
		}
                K += .5 * (velX * velX + velY * velY);
        }
}


number DifferentialAdhesion::f_interaction(BaseField *p, int q, BaseField *pp, int qq) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
        number dx = p->fieldDX[q];
        number dy = p->fieldDY[q];
	number xleft, xright, ybottom, ytop;
	
	//The part commented below serves for when the sublattice has dirichlet boundary conditions (e.g. zero at the boundaries)
	/*if(p->neighbors_sub[5+q*9]==-1)xright=0;
	else xright=p->fieldScalar[p->neighbors_sub[5+q*9]];

        if(p->neighbors_sub[7+q*9]==-1)ybottom=0;
        else ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]];

        if(p->neighbors_sub[3+q*9]==-1)xleft=0;
        else xleft=p->fieldScalar[p->neighbors_sub[3+q*9]];

        if(p->neighbors_sub[1+q*9]==-1)ytop=0;
        else ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];*/

	//this par gets the field values in teh respective directions from q;
        //It is hardcoded so take care, the relvant part is that the lattice is square;
        //The neighbors start form the top and rotate couterclockwise.
	xright=p->fieldScalar[p->neighbors_sub[5+q*9]]; 
	ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]]; 
	xleft=p->fieldScalar[p->neighbors_sub[3+q*9]]; 
	ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];

	number laplacianPhi = xright + ybottom + xleft + ytop - 4.*p->fieldScalar[q];

	// CH term coupled to chemical
	number CH =+ gamma*(8*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-2*p->fieldScalar[q])/lambda - 2*lambda*laplacianPhi);
   
	// area conservation term
	number A = - 4*mu/a0*(1-p->area/a0)*p->fieldScalar[q];

	// repulsion term
	number Rep = + 4*kappa/lambda*p->fieldScalar[q]*(phi2[k]-p->fieldScalar[q]*p->fieldScalar[q]);

	// adhesion term
	number suppress, laplacianSquare;
	number lsquare = 2 * p->fieldScalar[q] * laplacianPhi + 2 * (dx *dx + dy * dy);
	laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];
	suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));

	/*
	if(vec_omega[q + field_start_index[p->index]] < -0.5){
		laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];
		suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
	}
	else{
		if(phi_omega[q + field_start_index[p->index]] < -0.5)suppress = 0.;
		else {
			laplacianSquare = pp->fieldScalar[pp->neighbors_sub[5+qq*9]] + pp->fieldScalar[pp->neighbors_sub[7+qq*9]] + pp->fieldScalar[pp->neighbors_sub[3+qq*9]] +  pp->fieldScalar[pp->neighbors_sub[1+qq*9]] - 4.*pp->fieldScalar[qq];
			suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
		}
	}
	*/
	number Adh = - 4*lambda*omega*suppress*p->fieldScalar[q];

	// delta F / delta phi_i
	number V = CH + A + Rep + Adh;
	p->freeEnergy[q] += V;
	p->Pressure[q] = Rep - CH - A;

	return V;
}


void DifferentialAdhesion::calc_internal_forces(BaseField *p, int q) {

	int k = p->map_sub_to_box[q];

	//passive (passive force)
	number f_passive_x = 0;//(-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]] );
	number f_passive_y = 0;//(-1) * 0.5 * ( p->freeEnergy[p->neighbors_sub[7+q*9]] - p->freeEnergy[p->neighbors_sub[1+q*9]] );
	p->Fpassive_x[q] = f_passive_x;
	p->Fpassive_y[q] = f_passive_y;

	//active inter cells (active force)
	number fQ_self_x = -(p->S00*p->fieldDX[q] + p->S01*p->fieldDY[q]);
	number fQ_self_y = -(p->S01*p->fieldDX[q] - p->S00*p->fieldDY[q]);

	number fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ01[box->neighbors[7+k*9]] - sumQ01[box->neighbors[1+k*9]] ) ) - fQ_self_x;
	number fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( sumQ00[box->neighbors[7+k*9]] - sumQ00[box->neighbors[1+k*9]] ) ) - fQ_self_y;

	p->Factive_x[q] = zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive_y[q] = zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	p->velocityX[q] = (f_passive_x + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter) / friction;
	p->velocityY[q] = (f_passive_y + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter) / friction;
}


void DifferentialAdhesion::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	p->Q00 = p->S00;
	p->Q01 = p->S01;
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
	}

	//field_start_index[p->index]=size_rows;
	//size_rows += p->subSize;
}
