#include "ActiveNematic.h"
#include "../Utilities/Utils.h"

ActiveNematic::ActiveNematic() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				omega(0.004),
				mu(3.),
				kappa(0.1),
				friction(2.),
				zetaQ_self(0),
				zetaQ_inter(0),
				J_Q(1) {
	a0=PI*R*R;
}

ActiveNematic::~ActiveNematic() {

}


void ActiveNematic::get_settings(input_file &inp) {
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
}

void ActiveNematic::init() {
        a0=PI*R*R;
}

void ActiveNematic::read_topology(std::vector<BaseField*> &fields) {
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

void ActiveNematic::allocate_fields(std::vector<BaseField *> &fields) {
        for(int i = 0; i < (int) fields.size(); i++) {
                fields[i] = new MultiPhaseField();
        }
}

void ActiveNematic::apply_changes_after_equilibration(){
	zetaQ_self=zetaQ_self_active;
	zetaQ_inter=zetaQ_inter_active;
}

void ActiveNematic::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
	if(box->lees_edwards)throw RCexception("Interaction is not compatible with LEBc. Aborting");

	phi2.resize(Lx*Ly);
	sumQ00.resize(Lx*Ly);
	sumQ01.resize(Lx*Ly);
	for(int i =0; i<Lx*Ly; i++){resetSums(i);}
}

void ActiveNematic::resetSums(int k) {
	phi2[k]=0;
        sumQ00[k]=0;
        sumQ01[k]=0;
}


void ActiveNematic::updateFieldProperties(BaseField *p, int q, int k) {
	BaseInteraction::updateFieldProperties(p, q, k);
	number dx = p->fieldDX[q]; 
	number dy = p->fieldDY[q]; 
	p->S00 += -0.5*(dx*dx-dy*dy);
	p->S01 += -dx*dy;
}


void ActiveNematic::check_input_sanity(std::vector<BaseField *> &fields) {

}


void ActiveNematic::begin_energy_computation() {
		
        for(int i = 0; i < CONFIG_INFO->N(); i++) {
                initFieldProperties(CONFIG_INFO->fields()[i]);
        }
}

void ActiveNematic::initFieldProperties(BaseField *p) {

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


void ActiveNematic::begin_energy_computation(std::vector<BaseField *> &fields) {

	for(auto p : fields) {
		for(int q=0; q<p->subSize;q++)
			computeGlobalSums(p, q, false);
	}

        U = (number) 0;
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++)
			U += f_interaction(p, q);
        }

        K =0.;
	velX=0.;
	velY=0.;
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++){
			calc_internal_forces(p, q);
                	velX += p->Fpassive_x[q] + p->Factive_x[q];
                	velY += p->Fpassive_y[q] + p->Factive_y[q];
		}
                K += .5 * (velX * velX + velY * velY);
        }
	//std::cout<<U<<" "<<K<<std::endl;
	//exit(911);
}

void ActiveNematic::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;

	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
}

number ActiveNematic::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
        p->fieldDX[q] = dx;
        p->fieldDY[q] = dy;

	//std::cout<<k<<" "<<dx<<" "<<dy<<" "<<p->neighbors_sub[7+q*9]<<std::endl;

	//this part gets the field values in teh respective directions from q;
	//It is hardcoded so take care, the relvant part is that the lattice is square;
	//The neighbors start form the top and rotate couterclockwise.
	//number xleft, xright, ybottom, ytop;
	//xright=p->fieldScalar[p->neighbors_sub[5+q*9]]; 
	//ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]]; 
	//xleft=p->fieldScalar[p->neighbors_sub[3+q*9]]; 
	//ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];
	
	number laplacianPhi = p->fieldScalar[p->neighbors_sub[5+q*9]] + p->fieldScalar[p->neighbors_sub[7+q*9]] + p->fieldScalar[p->neighbors_sub[3+q*9]] + p->fieldScalar[p->neighbors_sub[1+q*9]] - 4.*p->fieldScalar[q];

	//xright=phi2[box->neighbors[5+k*9]]; 
	//ybottom=phi2[box->neighbors[7+k*9]]; 
	//xleft=phi2[box->neighbors[3+k*9]]; 
	//ytop=phi2[box->neighbors[1+k*9]];

	number laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];
	
	//if(p->index==0 && q==0)std::cout<<k<<" "<<phi2[k]<<" "<<phi2[box->neighbors[5+k*9]]<<" "<<phi2[box->neighbors[7+k*9]]<<" "<<phi2[box->neighbors[3+k*9]]<<" "<<phi2[box->neighbors[1+k*9]]  <<std::endl;

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

	//if(p->index==0 && q==0)std::cout<<CH<<" "<<A<<" "<<Rep<<" "<<Adh<<" "<< laplacianSquare<<" "<< lsquare <<std::endl;

	return V;
}


void ActiveNematic::calc_internal_forces(BaseField *p, int q) {

        //int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];

	//passive (passive force)
	p->Fpassive_x[q] = p->freeEnergy[q]*p->fieldDX[q];
	p->Fpassive_y[q] = p->freeEnergy[q]*p->fieldDY[q];

	//active inter cells (active force)
	number fQ_self_x = -(p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]);
	number fQ_self_y = -(p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]);

	number fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ01[box->neighbors[7+k*9]] - sumQ01[box->neighbors[1+k*9]] ) ) - fQ_self_x;
	number fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( sumQ00[box->neighbors[7+k*9]] - sumQ00[box->neighbors[1+k*9]] ) ) - fQ_self_y;

	p->Factive_x[q] = zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive_y[q] = zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	p->velocityX[q] = (p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter)/friction;
	p->velocityY[q] = (p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter)/friction;
}


void ActiveNematic::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	if(store)p->nemQ_old = {p->nemQ[0] , p->nemQ[1]};
	
	number t = 0.5 * atan2(p->S01, p->S00);
	std::vector<number> d = {cos(t) , sin(t)};
	number sgn = (d[0] * p->nemQ[0] + d[1] * p->nemQ[1] > 0.0)? 1.0:-1.0;
	p->nemQ[0] = p->nemQ_old[0] + dt * J_Q * (sgn * d[0] - p->nemQ[0]);
	p->nemQ[1] = p->nemQ_old[1] + dt * J_Q * (sgn * d[1] - p->nemQ[1]);

	p->Q00 = 0.5 * (p->nemQ[0] * p->nemQ[0] - p->nemQ[1] * p->nemQ[1]);
	p->Q01 = p->nemQ[0] * p->nemQ[1];

	if(anchoring) updateAnchoring(p);
}


void ActiveNematic::updateAnchoring(BaseField*p){

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
