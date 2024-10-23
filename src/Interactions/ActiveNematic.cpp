#include "ActiveNematic.h"
#include "../Utilities/Utils.h"

ActiveNematic::ActiveNematic() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				omega(0.004),
				mu(3.),
				R(8),
				kappa(0.1),
				friction(2.),
				zetaQ_self(0.002),
				zetaQ_inter(0.002),
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
	getInputNumber(&inp, "zetaQ_self", &zetaQ_self, 0);
	getInputNumber(&inp, "zetaQ_inter", &zetaQ_inter, 0);
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

void ActiveNematic::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
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


void ActiveNematic::updateFieldProperties(BaseField *p, int q) {
	BaseInteraction::updateFieldProperties(p, q);
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

	int sub=p->subSize;
	for(int q=0; q<sub;q++) {
		BaseInteraction::updateFieldProperties(p, q);
	        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
	        number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
	        p->fieldDX[q] = dx;
	        p->fieldDY[q] = dy;

		p->S00 += -0.5*(dx*dx-dy*dy);
		p->S01 += -dx*dy;
	}
}


void ActiveNematic::begin_energy_computation(std::vector<BaseField *> &fields) {

	for(auto p : fields) {
        	int sub=p->subSize;
		for(int q=0; q<sub;q++)
			computeGlobalSums(p, q, false);
	}

        U = (number) 0;
        for(auto p : fields) {
	        int sub=p->subSize;
                for(int q=0; q<sub;q++)
			U += f_interaction(p, q);
        }

        K = (number) 0;
        for(auto p : fields) {
        	int sub=p->subSize;
                p->Factive = std::vector<number> {0., 0.};
                p->Fpassive = std::vector<number> {0., 0.};
                for(int q=0; q<sub;q++)
			calc_internal_forces(p, q);
                velX = p->Fpassive[0] + p->Factive[0];
                velY = p->Fpassive[1] + p->Factive[1];
                K += .5 * (velX * velX + velY * velY);
        }

}

void ActiveNematic::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;
}

number ActiveNematic::f_interaction(BaseField *p, int q) {

	int  k  = p->GetSubIndex(q, box);
	number a = p->area;	
	number phi  = p->fieldScalar[q];

        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
        p->fieldDX[q] = dx;
        p->fieldDY[q] = dy;

	//this part gets the field values in teh respective directions from q;
	//It is hardcoded so take care, the relvant part is that the lattice is square;
	//The neighbors start form the top and rotate couterclockwise.
	number xleft, xright, ybottom, ytop;
	xright=p->fieldScalar[p->neighbors_sub[5+q*9]]; 
	ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]]; 
	xleft=p->fieldScalar[p->neighbors_sub[3+q*9]]; 
	ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];

	number laplacianPhi = xright + ybottom + xleft + ytop - 4.*p->fieldScalar[q];

	xright=phi2[box->neighbors[5+k*9]]; 
	ybottom=phi2[box->neighbors[7+k*9]]; 
	xleft=phi2[box->neighbors[3+k*9]]; 
	ytop=phi2[box->neighbors[1+k*9]];

	number laplacianSquare = xright + ybottom + xleft + ytop - 4.*phi2[k];

	// CH term coupled to chemical
	number CH =+ gamma*(8*phi*(1-phi)*(1-2*phi)/lambda - 2*lambda*laplacianPhi);
   
	// area conservation term
	number A = - 4*mu/a0*(1-a/a0)*phi;

	// repulsion term
	number Rep = + 4*kappa/lambda*phi*(phi2[k]-phi*phi);

	// adhesion term
	number lsquare = 2 * phi * laplacianPhi + 2 * (dx *dx + dy * dy);
	number suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
	number Adh = - 4*lambda*omega*suppress*phi;

	// delta F / delta phi_i
	number V = CH + A + Rep + Adh;
	p->freeEnergy[q] += V;

	return V;
}


void ActiveNematic::calc_internal_forces(BaseField *p, int q) {

        int  k  = p->GetSubIndex(q, box);
        // store derivatives
        //number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        //number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
        //p->fieldDX[q] = dx;
        //p->fieldDY[q] = dy;
        number dx = p->fieldDX[q];
        number dy = p->fieldDY[q];

	//passive (passive force)
	p->Fpassive[0] += p->freeEnergy[q]*dx;
	p->Fpassive[1] += p->freeEnergy[q]*dy;

	//active inter cells (active force)
	number fQ_self_x = -zetaQ_self*(p->Q00*dx + p->Q01*dy);
	number fQ_self_y = -zetaQ_self*(p->Q01*dx - p->Q00*dy);

	number fQ_inter_x = - zetaQ_inter * ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ01[box->neighbors[1+k*9]] - sumQ01[box->neighbors[7+k*9]] ) ) - fQ_self_x;
	number fQ_inter_y = - zetaQ_inter * ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) + 0.5 * ( sumQ00[box->neighbors[1+k*9]] - sumQ00[box->neighbors[7+k*9]] ) ) - fQ_self_y;

	p->Factive[0] += fQ_self_x + fQ_inter_x;
	p->Factive[1] += fQ_self_y + fQ_inter_y;

	p->velocityX[q] = (p->freeEnergy[q]*dx + fQ_self_x + fQ_inter_x)/friction;
	p->velocityY[q] = (p->freeEnergy[q]*dy + fQ_self_y + fQ_inter_y)/friction;
}


void ActiveNematic::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	/*if(store)
		p->thetaQ_old = p->thetaQ;

	number F00 = p->S00;
	number F01 = p->S01;

	number deformation = sqrt(sqrt(F01 * F01 + F00 * F00));

	p->thetaQ = p->thetaQ_old - dt * J_Q * deformation * atan2(F00 * p->Q01 - F01 * p->Q00, F00 * p->Q00 + F01 * p->Q01);
	p->Q00 = cos(2 * p->thetaQ);
	p->Q01 = sin(2 * p->thetaQ);*/

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

	number radius = R;
	number delta = -PI/12;
	number theta = 0;
	number S = 0.5;

	if(p->CoM[1] < 2.5 * radius){
		theta = (PI/2 + delta);
		p->Q00 = S * cos(2*theta);
		p->Q01 = S * sin(2*theta);
	}
	else if(p->CoM[1] > box->getYsize() - 2.5 * radius){
		theta = -(PI/2 + delta);
		p->Q00 = S * cos(2*theta);
		p->Q01 = S * sin(2*theta);
	}
}
