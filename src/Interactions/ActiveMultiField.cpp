#include "ActiveMultiField.h"
#include "../Utilities/Utils.h"

ActiveMultiField::ActiveMultiField() :
				BaseInteraction(),
				gamma(0.06),
				lambda(2.),
				mu(20.),
				kappa(0.5),
				friction(1.),
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
	if(box->lees_edwards)throw RCexception("Interaction is not compatible with LEBc. Aborting");
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
	p->S00 += -0.5*(p->fieldDX[q]*p->fieldDX[q]-p->fieldDY[q]*p->fieldDY[q]);
	p->S01 += -p->fieldDX[q]*p->fieldDY[q];
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
	        number dx = BaseInteraction::derivX(p, q, k);
	        number dy = BaseInteraction::derivY(p, q, k);
	        p->fieldDX[q] = dx;
	        p->fieldDY[q] = dy;

		p->S00 += -0.5*(dx*dx-dy*dy);
		p->S01 += -dx*dy;
	}
}


void ActiveMultiField::begin_energy_computation(std::vector<BaseField *> &fields) {

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
                	velX = p->Fpassive_x[q] + p->Factive_x[q];
                	velY = p->Fpassive_y[q] + p->Factive_y[q];
		}
                K += .5 * (velX * velX + velY * velY);
        }

}

void ActiveMultiField::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumS00[k]+=p->fieldScalar[q]*p->S00;
        sumS01[k]+=p->fieldScalar[q]*p->S01;
	BaseInteraction::update_sub_to_box_map(p, q, k, p->GetSubXIndex(q, box), p->GetSubYIndex(q, box));

        p->fieldDX[q] = BaseInteraction::derivX(p, q, k);
        p->fieldDY[q] = BaseInteraction::derivY(p, q, k);
	BaseInteraction::updateFieldProperties(p, q, k);
}

number ActiveMultiField::f_interaction(BaseField *p, int q) {

	int k = p->map_sub_to_box[q];
        //number dx = p->fieldDX[q];
        //number dy = p->fieldDY[q];
	//number xleft, xright, ybottom, ytop;
	
	//xright=p->fieldScalar[p->neighbors_sub[5+q*9]]; 
	//ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]]; 
	//xleft=p->fieldScalar[p->neighbors_sub[3+q*9]]; 
	//ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];

	number laplacianPhi = BaseInteraction::Laplacian(p, q, k);
	//xright + ybottom + xleft + ytop - 4.*p->fieldScalar[q];

	// CH term coupled to chemical
	number CH = gamma*(8*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-2*p->fieldScalar[q])/lambda - 2*lambda*laplacianPhi);
	number CH_density = (gamma/lambda)*4*p->fieldScalar[q]*p->fieldScalar[q]*(1-p->fieldScalar[q])*(1-p->fieldScalar[q]) + gamma*lambda*(p->fieldDX[q]*p->fieldDX[q] + p->fieldDY[q]*p->fieldDY[q]);
   
	// area conservation term
	number A = - 4*mu/a0*(1-p->area/a0)*p->fieldScalar[q];
	number A_density = - mu/a0*p->fieldScalar[q]*p->fieldScalar[q];

	// repulsion term
	number Rep = 4*kappa/lambda*p->fieldScalar[q]*(phi2[k]-p->fieldScalar[q]*p->fieldScalar[q]);
	number Rep_density = (kappa/lambda)*p->fieldScalar[q]*p->fieldScalar[q]*(phi2[k]-p->fieldScalar[q]*p->fieldScalar[q]);

	// adhesion term
	//number lsquare = 2 * p->fieldScalar[q] * laplacianPhi + 2 * (dx *dx + dy * dy);
	//number laplacianSquare = phi2[box->neighbors[5+k*9]] + phi2[box->neighbors[7+k*9]] + phi2[box->neighbors[3+k*9]] +  phi2[box->neighbors[1+k*9]] - 4.*phi2[k];
	//number suppress = (laplacianSquare-lsquare)/sqrt(1+(laplacianSquare-lsquare)*(laplacianSquare-lsquare));
	//number Adh = - 4*lambda*omega*suppress*p->fieldScalar[q];
	number Adh = 0.;

	// delta F / delta phi_i
	number V = CH + A + Rep + Adh;
	p->freeEnergy[q] += V;
	p->freeEnergyDensity[q] += CH_density + A_density + Rep_density;
	p->freeEnergyDensityGradient_x[q] = 2 * gamma * lambda * p->fieldDX[q];
	p->freeEnergyDensityGradient_y[q] = 2 * gamma * lambda * p->fieldDY[q];
	p->Pressure[q] = Rep - CH - A;

	return V;
}


void ActiveMultiField::calc_internal_forces(BaseField *p, int q) {

        int  k  = p->GetSubIndex(q, box);
	//passive (passive force)
	p->Fpassive_x[q] = p->freeEnergy[q]*p->fieldDX[q];
	p->Fpassive_y[q] = p->freeEnergy[q]*p->fieldDY[q];

	//active inter cells (active force)
	p->Factive_x[q] = zetaS*sumS00[k]*p->fieldDX[q] + zetaS*sumS01[k]*p->fieldDY[q];
	p->Factive_y[q] = zetaS*sumS01[k]*p->fieldDX[q] - zetaS*sumS00[k]*p->fieldDY[q];


	if(box->getWalls(k)<wall_slip){
		p->velocityX[q] = (p->freeEnergy[q]*p->fieldDX[q] - zetaS*sumS00[k]*p->fieldDX[q] - zetaS*sumS01[k]*p->fieldDY[q])/friction;
		p->velocityY[q] = (p->freeEnergy[q]*p->fieldDY[q] - zetaS*sumS01[k]*p->fieldDX[q] + zetaS*sumS00[k]*p->fieldDY[q])/friction;
	}
	else{
		p->velocityX[q] = 0.;
		p->velocityY[q] = 0.;
	}
}
