#include "LEBcActiveNematic.h"
#include "../Utilities/Utils.h"

LEBcActiveNematic::LEBcActiveNematic() :
				BaseInteraction(),
				gamma(0.01),
				lambda(2.5),
				omega(0.004),
				mu(3.),
				kappa(0.1),
				friction(2.),
				zetaQ_self(0),
				zetaQ_inter(0),
				J_Q(1),
				shear_rate(0.) {
	a0=PI*R*R;
}

LEBcActiveNematic::~LEBcActiveNematic() {

}


void LEBcActiveNematic::get_settings(input_file &inp) {
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
	getInputNumber(&inp, "lees_edwards_shear_rate", &shear_rate_active, 1);
}

void LEBcActiveNematic::init() {
        a0=PI*R*R;
}

void LEBcActiveNematic::read_topology(std::vector<BaseField*> &fields) {
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

void LEBcActiveNematic::allocate_fields(std::vector<BaseField *> &fields) {
        for(int i = 0; i < (int) fields.size(); i++) {
                fields[i] = new LEBcMultiPhaseField();
        }
}

void LEBcActiveNematic::apply_changes_after_equilibration(){
	zetaQ_self=zetaQ_self_active;
	zetaQ_inter=zetaQ_inter_active;
	shear_rate = shear_rate_active;
}

void LEBcActiveNematic::set_box(BaseBox *boxArg) {
	box = boxArg;
	int Lx=box->getXsize();
	int Ly=box->getYsize();
	if(!box->lees_edwards)throw RCexception("LEBc interaction...but no LEBc box...Aborting");

	phi2.resize(Lx*Ly);
	sumQ00.resize(Lx*Ly);
	sumQ01.resize(Lx*Ly);
	for(int i =0; i<Lx*Ly; i++){resetSums(i);}
}

void LEBcActiveNematic::resetSums(int k) {
	phi2[k]=0;
        sumQ00[k]=0;
        sumQ01[k]=0;
}


void LEBcActiveNematic::updateFieldProperties(BaseField *p, int q, int k) {
	BaseInteraction::updateFieldProperties(p, q, k);
	number dx = p->fieldDX[q]; 
	number dy = p->fieldDY[q]; 
	p->S00 += -0.5*(dx*dx-dy*dy);
	p->S01 += -dx*dy;
}


void LEBcActiveNematic::check_input_sanity(std::vector<BaseField *> &fields) {

}


void LEBcActiveNematic::begin_energy_computation() {
		
        for(int i = 0; i < CONFIG_INFO->N(); i++) {
                initFieldProperties(CONFIG_INFO->fields()[i]);
        }
}

void LEBcActiveNematic::initFieldProperties(BaseField *p) {

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


void LEBcActiveNematic::begin_energy_computation(std::vector<BaseField *> &fields) {

	box->setNeighborsPeriodic(box->getXsize(), box->getYsize());
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
        for(auto p : fields) {
                p->Factive = std::vector<number> {0., 0.};
                p->Fpassive = std::vector<number> {0., 0.};
                for(int q=0; q<p->subSize;q++)
			calc_internal_forces(p, q);
                velX = p->Fpassive[0] + p->Factive[0];
                velY = p->Fpassive[1] + p->Factive[1];
                K += .5 * (velX * velX + velY * velY);
        }

}

void LEBcActiveNematic::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	int k = p->GetSubIndex(q, box);
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;

	//if(q==0 || q==25 || q==26)std::cout<<"interaction: "<<k<<" "<<  box->getElementX(k, 0) <<" "<< box->getElementY(k, 0) <<std::endl;
	BaseInteraction::update_sub_to_box_map(p, q, k, box->getElementX(k, 0), box->getElementY(k, 0));
}

number LEBcActiveNematic::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
	number ytop=(box->weight_site[7+k*9]*p->fieldScalar[p->neighbors_sub[7+q*9]]+box->weight_site_next[7+k*9]*p->fieldScalar[p->neighbors_sub[6+q*9]]);
	number ybottom=(box->weight_site[1+k*9]*p->fieldScalar[p->neighbors_sub[1+q*9]]+box->weight_site_next[1+k*9]*p->fieldScalar[p->neighbors_sub[2+q*9]]);

        number dx = .5*( p->fieldScalar[p->neighbors_sub[5+q*9]] - p->fieldScalar[p->neighbors_sub[3+q*9]] );
        number dy = .5*( ytop - ybottom ); 

        p->fieldDX[q] = dx;
        p->fieldDY[q] = dy;

	//this part gets the field values in teh respective directions from q;
	//It is hardcoded so take care, the relvant part is that the lattice is square;
	//The neighbors start form the bottom and rotate couterclockwise.
	//number xleft, xright, ybottom, ytop;
	//xright=p->fieldScalar[p->neighbors_sub[5+q*9]]; 
	//ybottom=p->fieldScalar[p->neighbors_sub[7+q*9]]; 
	//xleft=p->fieldScalar[p->neighbors_sub[3+q*9]]; 
	//ytop=p->fieldScalar[p->neighbors_sub[1+q*9]];
	
	number laplacianPhi = p->fieldScalar[p->neighbors_sub[5+q*9]] + ytop + p->fieldScalar[p->neighbors_sub[3+q*9]] + ybottom - 4.*p->fieldScalar[q];

	//xright=phi2[box->neighbors[5+k*9]]; 
	//ybottom=phi2[box->neighbors[7+k*9]]; 
	//xleft=phi2[box->neighbors[3+k*9]]; 
	//ytop=phi2[box->neighbors[1+k*9]];
	
	ytop=(box->weight_site[7+k*9]*phi2[box->neighbors[7+k*9]]+box->weight_site_next[7+k*9]*phi2[box->neighbors[6+k*9]]);
	ybottom=(box->weight_site[1+k*9]*phi2[box->neighbors[1+k*9]]+box->weight_site_next[1+k*9]*phi2[box->neighbors[2+k*9]]);
	
	//if(p->index==0 && q==0)std::cout<< box->weight_site[7+k*9] <<" "<< phi2[box->neighbors[7+q*9]] <<" "<< box->weight_site_next[7+k*9] <<" "<< phi2[box->neighbors[8+q*9]] <<" "<<box->weight_site[1+k*9]  <<" "<< box->weight_site_next[1+k*9] << " "<<k <<" "<<phi2[k]<< std::endl;
	//if(p->index==0 && q==0)std::cout<<k<<" "<<phi2[k]<<" "<<phi2[box->neighbors[5+k*9]]<<" "<<phi2[box->neighbors[7+k*9]]<<" "<<phi2[box->neighbors[3+k*9]]<<" "<<phi2[box->neighbors[1+k*9]]  <<" "<<ytop<<" "<<ybottom<<" "<< box->weight_site[1+k*9] <<" " << box->weight_site_next[1+k*9] <<std::endl;

	number laplacianSquare = phi2[box->neighbors[5+k*9]] + ytop + phi2[box->neighbors[3+k*9]] +  ybottom - 4.*phi2[k];

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

	//if(p->index==0 && q==0)std::cout<<CH<<" "<<A<<" "<<Rep<<" "<<Adh<<" "<< laplacianSquare<<" "<< lsquare <<std::endl;

	return V;
}


void LEBcActiveNematic::calc_internal_forces(BaseField *p, int q) {

        //int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];

	//passive (passive force)
	p->Fpassive[0] += p->freeEnergy[q]*p->fieldDX[q];
	p->Fpassive[1] += p->freeEnergy[q]*p->fieldDY[q];

	//active inter cells (active force)
	number fQ_self_x = -(p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]);
	number fQ_self_y = -(p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]);

	number fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( (box->weight_site[7+k*9]*sumQ01[box->neighbors[7+q*9]]+box->weight_site_next[7+k*9]*sumQ01[box->neighbors[6+q*9]]) - (box->weight_site[1+k*9]*sumQ01[box->neighbors[1+q*9]]+box->weight_site_next[1+k*9]*sumQ01[box->neighbors[2+q*9]]) ) ) - fQ_self_x;

	number fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( (box->weight_site[7+k*9]*sumQ00[box->neighbors[7+q*9]]+box->weight_site_next[7+k*9]*sumQ00[box->neighbors[6+q*9]]) - (box->weight_site[1+k*9]*sumQ00[box->neighbors[1+q*9]]+box->weight_site_next[1+k*9]*sumQ00[box->neighbors[2+q*9]]) ) ) - fQ_self_y;

	p->Factive[0] += zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive[1] += zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	p->velocityX[q] = (p->freeEnergy[q]*p->fieldDX[q] + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter)/friction;
	p->velocityY[q] = (p->freeEnergy[q]*p->fieldDY[q] + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter)/friction; //- 0.05;
}


void LEBcActiveNematic::updateDirectedActiveForces(number dt, BaseField*p, bool store){

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


void LEBcActiveNematic::updateAnchoring(BaseField*p){

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

number LEBcActiveNematic::get_velocity_x(BaseField *p, int q){
	int y = p->map_sub_to_box_y[q];
	int y_sub = p->sub_corner_bottom_left / box->getXsize();
	int unrap_number = p->unrap_sub_corner_bottom_left_y / box->getYsize();

	if(p->unrap_sub_corner_bottom_left_y>=0){
		//if(q==p->subSize-1)std::cout<<"velocity: "<<unrap_number<<" "<< p->velocityX[q] << " " << unrap_number * shear_rate * box->getYsize() << " "<< p->velocityX[q] - unrap_number * shear_rate * box->getYsize() << " "<< y << " "<< box->getYsize()-y_sub  << std::endl;
		if(y<box->getYsize() && y>=y_sub)return p->velocityX[q] - unrap_number * shear_rate * box->getYsize();
		else return p->velocityX[q] - (unrap_number+1) * shear_rate * box->getYsize();
	}
	else{
		unrap_number-=1;
		//if(q==p->subSize-1)std::cout<<"velocity: "<<unrap_number<<" "<< p->velocityX[q] << " " << (number)unrap_number * shear_rate * box->getYsize() << " "<< p->velocityX[q] - unrap_number * shear_rate * box->getYsize() << " "<< y << " "<< box->getYsize() <<" "<<y_sub << std::endl;
		if(y<box->getYsize() && y>=y_sub)return p->velocityX[q] - unrap_number * shear_rate * box->getYsize();
		else return p->velocityX[q] - (unrap_number+1) * shear_rate * box->getYsize();
	}
}

number LEBcActiveNematic::get_velocity_y(BaseField *p, int q){return p->velocityY[q];}
