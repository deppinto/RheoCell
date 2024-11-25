#include "SimpleMultiField.h"
#include "../Utilities/Utils.h"

SimpleMultiField::SimpleMultiField() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				mu(3.),
				kappa(0.1) {
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

	number V;
	//int Lx=p->LsubX;
        //int Ly=p->LsubY;
	number a = p->area;	
	int  k  = p->GetSubIndex(q, box);
	//std::vector<int> s = p->neighbors_sub;
	//std::vector<number> f = p->fieldScalar;
	number phi  = p->fieldScalar[q];
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

	// delta F / delta phi_i
	V = (

	// CH term coupled to chemical
	+ gamma*(8*phi*(1-phi)*(1-2*phi)/lambda - 2*lambda*laplacianPhi)
   
	// area conservation term
	- 4*mu/a0*(1-a/a0)*phi

	// repulsion term
	+ 4*kappa/lambda*phi*(phi2[k]-phi*phi)
	);

	p->freeEnergy[q] = V;
	return V;
}


void SimpleMultiField::begin_energy_computation() {

}
