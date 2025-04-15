#include "SimpleMultiField.h"
#include "../Utilities/Utils.h"

SimpleMultiField::SimpleMultiField() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				mu(3.),
				kappa(0.1),
				omega(0.) {
	a0=PI*R*R;
}

SimpleMultiField::~SimpleMultiField() {

}


void SimpleMultiField::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);

	getInputInt(&inp, "R", &R, 0);
	getInputNumber(&inp, "gamma", &gamma, 0);
	getInputNumber(&inp, "mu", &mu, 0);
	getInputNumber(&inp, "kappa", &kappa, 0);
	getInputNumber(&inp, "omega", &omega, 0);
}


void SimpleMultiField::read_topology(std::vector<BaseField *> &fields) {
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


void SimpleMultiField::set_box(BaseBox *boxArg) {
                box = boxArg;
		int Lx=box->getXsize();
	        int Ly=box->getYsize();
		if(box->lees_edwards)throw RCexception("Interaction is not compatible with LEBc. Aborting");
        	phi2.resize(Lx*Ly);
}

void SimpleMultiField::init() {
	a0=PI*R*R;
}

void SimpleMultiField::allocate_fields(std::vector<BaseField *> &fields) {
	for(int i = 0; i < (int) fields.size(); i++) {
		fields[i] = new MultiPhaseField();
	}
}

void SimpleMultiField::check_input_sanity(std::vector<BaseField *> &fields) {

}

void SimpleMultiField::resetSums(int k) {
	phi2[k]=0;
}

void SimpleMultiField::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));
}

void SimpleMultiField::begin_energy_computation(std::vector<BaseField *> &fields) {

	box->UpdateWalls(false);
	for(auto p : fields) {
		for(int q=0; q<p->subSize;q++)
			computeGlobalSums(p, q, false);
	}

        U = (number) 0;
        for(auto p : fields) {
                for(int q=0; q<p->subSize;q++) {
                        U += f_interaction(p, q);
                }
        }

        K = (number) 0;
}

number SimpleMultiField::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( p->fieldScalar[p->neighbors_sub[7+q*9]] - p->fieldScalar[p->neighbors_sub[1+q*9]] );
        p->fieldDX[q] = dx;
        p->fieldDY[q] = dy;
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
	p->Pressure[q] = Rep - CH - A;

	return V;
}


void SimpleMultiField::begin_energy_computation() {

}
