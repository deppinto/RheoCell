#include "ActiveMultiField.h"
#include "../Utilities/Utils.h"

ActiveMultiField::ActiveMultiField() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				mu(3.),
				R(8),
				kappa(0.1),
				friction(2.),
				zetaS(0) {
	a0=PI*R*R;
}

ActiveMultiField::~ActiveMultiField() {

}


void ActiveMultiField::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "R", &R, 0);
	getInputNumber(&inp, "lambda", &lambda, 0);
	getInputNumber(&inp, "gamma", &gamma, 0);
	getInputNumber(&inp, "mu", &mu, 0);
	getInputNumber(&inp, "kappa", &kappa, 0);
	getInputNumber(&inp, "friction", &friction, 0);
	getInputNumber(&inp, "zetaS", &zetaS_active, 0);
}

void ActiveMultiField::init() {
        a0=PI*R*R;
}

void ActiveMultiField::read_topology(std::vector<BaseField*> &fields) {
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

void ActiveMultiField::allocate_fields(std::vector<BaseField *> &fields) {
        for(int i = 0; i < (int) fields.size(); i++) {
                fields[i] = new MultiPhaseField();
        }
}

void ActiveMultiField::apply_changes_after_equilibration(){
	zetaS=zetaS_active;
}

void ActiveMultiField::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
	phi2.resize(Lx*Ly);
	sumS00.resize(Lx*Ly);
	sumS01.resize(Lx*Ly);
	for(int i =0; i<Lx*Ly; i++){resetSums(i);}
}

void ActiveMultiField::resetSums(int k) {
	phi2[k]=0;
        sumS00[k]=0;
        sumS01[k]=0;
}


void ActiveMultiField::updateFieldProperties(BaseField *p, int q, int k) {
	BaseInteraction::updateFieldProperties(p, q, k);
	number dx = p->fieldDX[q]; 
	number dy = p->fieldDY[q]; 
	p->S00 += -0.5*(dx*dx-dy*dy);
	p->S01 += -dx*dy;
}


void ActiveMultiField::check_input_sanity(std::vector<BaseField *> &fields) {

}


void ActiveMultiField::begin_energy_computation() {
		
        for(int i = 0; i < CONFIG_INFO->N(); i++) {
                initFieldProperties(CONFIG_INFO->fields()[i]);
        }
}

void ActiveMultiField::initFieldProperties(BaseField *p) {

	int sub=p->subSize;
	for(int q=0; q<sub;q++) {
		int k = p->GetSubIndex(q, box);
		BaseInteraction::updateFieldProperties(p, q, k);
	        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
	        number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
	        p->fieldDX[q] = dx;
	        p->fieldDY[q] = dy;

		p->S00 += -0.5*(dx*dx-dy*dy);
		p->S01 += -dx*dy;
	}
}


void ActiveMultiField::begin_energy_computation(std::vector<BaseField *> &fields) {

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

void ActiveMultiField::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumS00[k]+=p->fieldScalar[q]*p->S00;
        sumS01[k]+=p->fieldScalar[q]*p->S01;
	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
}

number ActiveMultiField::f_interaction(BaseField *p, int q) {

	int  k  = p->GetSubIndex(q, box);
	number a = p->area;	
	number phi  = p->fieldScalar[q];

	//this part gets the field values in teh respective directions from q;
	//It is hardcoded so take care, the relvant part is that the lattice is square;
	//The neighbors start form the top and rotate couterclockwise.
	number xleft, xright, ybottom, ytop;
	xright=p->fieldScalar[p->neighbors_sub[5+q*9]]; 
	ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]]; 
	xleft=p->fieldScalar[p->neighbors_sub[3+q*9]]; 
	ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];

	number laplacianPhi = xright + ybottom + xleft + ytop - 4.*p->fieldScalar[q];

	// CH term coupled to chemical
	number CH =+ gamma*(8*phi*(1-phi)*(1-2*phi)/lambda - 2*lambda*laplacianPhi);
   
	// area conservation term
	number A = - 4*mu/a0*(1-a/a0)*phi;

	// repulsion term
	number Rep = + 4*kappa/lambda*phi*(phi2[k]-phi*phi);

	// delta F / delta phi_i
	number V = CH + A + Rep;
	p->freeEnergy[q] += V;

	return V;
}


void ActiveMultiField::calc_internal_forces(BaseField *p, int q) {

        int  k  = p->GetSubIndex(q, box);
        // store derivatives
        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( p->fieldScalar[p->neighbors_sub[1+q*9]] - p->fieldScalar[p->neighbors_sub[7+q*9]] );
        p->fieldDX[q] = dx;
        p->fieldDY[q] = dy;

	//passive (passive force)
	p->Fpassive[0] += p->freeEnergy[q]*dx;
	p->Fpassive[1] += p->freeEnergy[q]*dy;

	//active inter cells (active force)
	p->Factive[0] += zetaS*sumS00[k]*dx + zetaS*sumS01[k]*dy;
	p->Factive[1] += zetaS*sumS01[k]*dx - zetaS*sumS00[k]*dy;

	p->velocityX[q] = (p->freeEnergy[q]*dx - zetaS*sumS00[k]*dx - zetaS*sumS01[k]*dy)/friction;
	p->velocityY[q] = (p->freeEnergy[q]*dy - zetaS*sumS01[k]*dx + zetaS*sumS00[k]*dy)/friction;
}
