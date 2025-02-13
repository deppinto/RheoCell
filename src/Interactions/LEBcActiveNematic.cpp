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
	p->shear_velocity_sign[q] = 0;
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
		//int k = p->GetSubIndex(q, box);
		int k = p->map_sub_to_box[q];
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

		p->velocityX_CoM = 0.;
		p->phi_CoM = 0.;
		std::fill(p->velocityX_correction.begin(), p->velocityX_correction.end(), 0);
		std::fill(p->phi_correction.begin(), p->phi_correction.end(), 0);
		if(p->LsubY + int(p->sub_corner_bottom_left/box->getXsize()) >= box->getYsize())p->update_positions(box);

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
                for(int q=0; q<p->subSize;q++){
			calc_internal_forces(p, q);
                	velX = p->Fpassive_x[q] + p->Factive_x[q];
                	velY = p->Fpassive_y[q] + p->Factive_y[q];
		}
                K += .5 * (velX * velX + velY * velY);
        }

}

void LEBcActiveNematic::computeGlobalSums(BaseField *p, int q, bool update_global_sums) {

	//int k = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
	phi2[k]+=p->fieldScalar[q]*p->fieldScalar[q];
	sumQ00[k]+=p->fieldScalar[q]*p->Q00;
        sumQ01[k]+=p->fieldScalar[q]*p->Q01;

	//if(q==0 || q==25 || q==26)std::cout<<"interaction: "<<k<<" "<<  box->getElementX(k, 0) <<" "<< box->getElementY(k, 0) <<std::endl;
	//BaseInteraction::update_sub_to_box_map(p, q, k, box->getElementX(k, 0), box->getElementY(k, 0));
}

number LEBcActiveNematic::f_interaction(BaseField *p, int q) {

	//int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];
	//int ky = p->map_sub_to_box_y[q];

	number ytop, ybottom;
	/*if(ky==box->getYsize()-1){

		int qx = (int(q - int(q/p->LsubX) * p->LsubX) - int(box->get_shear_displacement()));
		int qy = (int(q/p->LsubX) + 1);
		while(qx<0)qx+=p->LsubX;
		if(qy>=p->LsubY)qy-=p->LsubY;
		//int delta_q = ((int(q - int(q/p->LsubX) * p->LsubX) + dist)%p->LsubX) + ((int(q/p->LsubX) + 1)%p->LsubY) * p->LsubX;
		int delta_q = qx + qy * p->LsubX;
		ytop=(box->weight_site[7+k*9]*p->fieldScalar[delta_q]+box->weight_site_next[7+k*9]*p->fieldScalar[p->neighbors_sub[3+delta_q*9]]);
		//ytop=(0.5*p->fieldScalar[delta_q]+0.5*p->fieldScalar[p->neighbors_sub[3+delta_q*9]]);
		ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
		
		if(p->index==0 && k - ky * box->getXsize() == 58)std::cout<<"Free Energy-----TOP: "<< k<<" "<<box->weight_site[7+k*9] << " "<<box->neighbors[7+k*9] - int(box->neighbors[7+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->weight_site_next[7+k*9]<<" "<<box->neighbors[6+k*9] - int(box->neighbors[6+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->neighbors[7+k*9]/box->getXsize()<<" "<<box->get_shear_displacement()  <<" "<<delta_q << " "<<p->fieldScalar[delta_q]<<" "<< p->fieldScalar[p->neighbors_sub[3+delta_q*9]]<<" "<<ytop<<" "<<ybottom<<" "<<p->map_sub_to_box[delta_q]<<" "<<p->map_sub_to_box[p->neighbors_sub[3+delta_q*9]]  <<std::endl;

	}
	else if(ky==0){

		int qx = (int(q - int(q/p->LsubX) * p->LsubX) + int(box->get_shear_displacement()));
		int qy = (int(q/p->LsubX) - 1);
		while(qx>=p->LsubX)qx-=p->LsubX;
		if(qy<0)qy+=p->LsubY;
		//int delta_q = ((int(q - int(q/p->LsubX) * p->LsubX) + dist)%p->LsubX) + ((int(q/p->LsubX) - 1)%p->LsubY) * p->LsubX;
		int delta_q = qx + qy * p->LsubX;
		ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
		ybottom=(box->weight_site[1+k*9]*p->fieldScalar[delta_q]+box->weight_site_next[1+k*9]*p->fieldScalar[p->neighbors_sub[5+delta_q*9]]);
		//ybottom=(0.5*p->fieldScalar[delta_q]+0.5*p->fieldScalar[p->neighbors_sub[5+delta_q*9]]);

		//if(p->index==97 && k - ky * box->getXsize() == 75)std::cout<<"Free Energy-----BOTTOM: "<< k<< " "<< box->weight_site[1+k*9] << " "<<box->neighbors[1+k*9] - int(box->neighbors[1+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->weight_site_next[1+k*9]<<" "<<box->neighbors[3+k*9] - int(box->neighbors[3+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->neighbors[1+k*9]/box->getXsize()<<" "<<box->get_shear_displacement()  <<std::endl;
	}
	else{
		ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
		ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
	}*/

	ytop=(box->weight_site[7+k*9]*p->fieldScalar[p->neighbors_sub[7+q*9]]+box->weight_site_next[7+k*9]*p->fieldScalar[p->neighbors_sub[6+q*9]]);
	ybottom=(box->weight_site[1+k*9]*p->fieldScalar[p->neighbors_sub[1+q*9]]+box->weight_site_next[1+k*9]*p->fieldScalar[p->neighbors_sub[2+q*9]]);
	
	//ytop=p->fieldScalar[p->neighbors_sub[7+q*9]];
	//ybottom=p->fieldScalar[p->neighbors_sub[1+q*9]];
	
	//if(ky==box->getYsize()-1 && p->index==0 && k - ky * box->getXsize() == 58)std::cout<<"Free Energy-----TOP: "<<q<<" "<< k<<" "<<box->weight_site[7+k*9] << " "<<box->neighbors[7+k*9] - int(box->neighbors[7+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->weight_site_next[7+k*9]<<" "<<box->neighbors[6+k*9] - int(box->neighbors[6+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->neighbors[7+k*9]/box->getXsize()<<" "<<box->get_shear_displacement()  <<" "<< p->neighbors_sub[7+q*9] <<" "<< p->neighbors_sub[6+q*9] <<" "<< p->neighbors_sub[1+q*9] <<" "<< p->neighbors_sub[2+q*9]  <<" "<<ytop<<" "<<ybottom<<" "<<p->map_sub_to_box[p->neighbors_sub[7+q*9]]<<" "<<p->map_sub_to_box[p->neighbors_sub[6+q*9]]  <<std::endl;


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

	//if(q==509 && p->index==1)std::cout<<"Free Energy 509--: "<< p->fieldScalar[p->neighbors_sub[5+q*9]]<<" "<< ytop<<" "<< p->fieldScalar[p->neighbors_sub[3+q*9]]<<" "<< ybottom <<" "<< 4.*p->fieldScalar[q] << std::endl;
	//if(ky==99 && p->index==0 && k - ky * box->getXsize() == 58)std::cout<<"Free Energy-----TOP: "<< box->weight_site[1+k*9] << " "<<box->neighbors[1+k*9] - int(box->neighbors[1+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->weight_site_next[1+k*9]<<" "<<box->neighbors[2+k*9] - int(box->neighbors[2+k*9]/box->getXsize()) * box->getXsize()<<" "<<box->neighbors[7+k*9]/box->getXsize()<<" "<<box->get_shear_displacement()  <<std::endl;

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

	//if(q==p->subSize-1 && p->index==1)std::cout<<"Free Energy S: "<<CH<<" "<<A<<" "<<Rep<<" "<<Adh<<" "<< laplacianSquare<<" "<< lsquare <<" "<<p->fieldScalar[q]<<" "<<dx<<" "<<dy<<" "<<k<<" "<<laplacianPhi<< std::endl;
	//if(p->map_sub_to_box_y[q]==99 && p->map_sub_to_box_x[q]==58  && p->index==0)std::cout<<"Free Energy 0: "<<q<<" "<<CH<<" "<<A<<" "<<Rep<<" "<<Adh<<" "<< laplacianSquare<<" "<< lsquare <<" "<<p->fieldScalar[q]<<" "<<dx<<" "<<dy<<" "<<k<<" "<<laplacianPhi<< std::endl;
	//if(q==509 && p->index==1)std::cout<<"Free Energy 509: "<<CH<<" "<<A<<" "<<Rep<<" "<<Adh<<" "<< laplacianSquare<<" "<< lsquare <<" "<<p->fieldScalar[q]<<" "<<dx<<" "<<dy<<" "<<k<<" "<<laplacianPhi<<" "<<ytop<<" " <<ybottom << std::endl;

	return V;
}


void LEBcActiveNematic::calc_internal_forces(BaseField *p, int q) {

        //int  k  = p->GetSubIndex(q, box);
	int k = p->map_sub_to_box[q];

	//passive (passive force)
	number f_passive_x = 0.;
	number f_passive_y = 0.;
	//f_passive_x = p->freeEnergy[q]*p->fieldDX[q];
	//f_passive_y = p->freeEnergy[q]*p->fieldDY[q];
	f_passive_x = (-1) * 0.5 * (p->freeEnergy[p->neighbors_sub[5+q*9]] - p->freeEnergy[p->neighbors_sub[3+q*9]]);
	f_passive_y = (-1) * 0.5 * ( (box->weight_site[7+k*9]*p->freeEnergy[p->neighbors_sub[7+q*9]]+box->weight_site_next[7+k*9]*p->freeEnergy[p->neighbors_sub[6+q*9]]) - (box->weight_site[1+k*9]*p->freeEnergy[p->neighbors_sub[1+q*9]]+box->weight_site_next[1+k*9]*p->freeEnergy[p->neighbors_sub[2+q*9]]) );

	p->Fpassive_x[q] = f_passive_x;
	p->Fpassive_y[q] = f_passive_y;

	//active inter cells (active force)
	number fQ_self_x = 0.;
	number fQ_self_y = 0.;

	if(zetaQ_self!=0){
		fQ_self_x = -(p->Q00*p->fieldDX[q] + p->Q01*p->fieldDY[q]);
		fQ_self_y = -(p->Q01*p->fieldDX[q] - p->Q00*p->fieldDY[q]);
	}

	number fQ_inter_x = 0.;
	number fQ_inter_y = 0.;

	if(zetaQ_inter!=0){
		fQ_inter_x = - ( 0.5 * ( sumQ00[box->neighbors[5+k*9]] - sumQ00[box->neighbors[3+k*9]] ) + 0.5 * ( (box->weight_site[7+k*9]*sumQ01[box->neighbors[7+q*9]]+box->weight_site_next[7+k*9]*sumQ01[box->neighbors[6+q*9]]) - (box->weight_site[1+k*9]*sumQ01[box->neighbors[1+q*9]]+box->weight_site_next[1+k*9]*sumQ01[box->neighbors[2+q*9]]) ) ) - fQ_self_x;

		fQ_inter_y = - ( 0.5 * ( sumQ01[box->neighbors[5+k*9]] - sumQ01[box->neighbors[3+k*9]] ) - 0.5 * ( (box->weight_site[7+k*9]*sumQ00[box->neighbors[7+q*9]]+box->weight_site_next[7+k*9]*sumQ00[box->neighbors[6+q*9]]) - (box->weight_site[1+k*9]*sumQ00[box->neighbors[1+q*9]]+box->weight_site_next[1+k*9]*sumQ00[box->neighbors[2+q*9]]) ) ) - fQ_self_y;
	}

	p->Factive_x[q] = zetaQ_self * fQ_self_x + zetaQ_inter * fQ_inter_x;
	p->Factive_y[q] = zetaQ_self * fQ_self_y + zetaQ_inter * fQ_inter_y;

	p->velocityX[q] = (f_passive_x + fQ_self_x * zetaQ_self + fQ_inter_x * zetaQ_inter)/friction;
	p->velocityY[q] = (f_passive_y + fQ_self_y * zetaQ_self + fQ_inter_y * zetaQ_inter)/friction; // - 0.1; 

	p->velocityX_correction[int(q/p->LsubX)] += p->velocityX[q] + shear_rate * (((double)p->map_sub_to_box_y[q] + 0.5) - 0.5 * (double)box->getYsize());
	p->phi_correction[int(q/p->LsubX)] += p->fieldScalar[q];
	p->velocityX_CoM += p->velocityX[q] + shear_rate * (((double)p->map_sub_to_box_y[q] + 0.5) - 0.5 * (double)box->getYsize());
	p->phi_CoM += p->fieldScalar[q];
;
}


void LEBcActiveNematic::updateDirectedActiveForces(number dt, BaseField*p, bool store){

	if(J_Q==0)return;

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

	//if(p->index==0 && p->map_sub_to_box_x[q]==50 && p->map_sub_to_box_y[q]==99)std::cout<<"average velocity along line: "<< p->map_sub_to_box_y[q]<<" "<< p->velocityX_correction[int(q/p->LsubX)] / p->phi_correction[int(q/p->LsubX)] <<" "<<shear_rate * (((double)p->map_sub_to_box_y[q] + 0.5) - 0.5 * (double)box->getYsize())<<" "<< p->velocityX[q]  <<" "<< p->velocityX_CoM / p->phi_CoM<<" "<<p->freeEnergy[q]<<" "<<p->fieldDX[q]<<" "<<p->fieldDY[q] <<std::endl;

	//return p->velocityX[q] + p->shear_velocity_sign[q] * shear_rate * box->getYsize();
	return p->velocityX[q] + shear_rate * (((double)p->map_sub_to_box_y[q] + 0.5) - 0.5 * (double)box->getYsize());
	//return p->velocityX[q] + p->shear_velocity_sign[q] * shear_rate * (box->getYsize() - 1) + shear_rate * (((double)p->map_sub_to_box_y[q] + 0.5) - 0.5 * (double)box->getYsize());
}

number LEBcActiveNematic::get_velocity_y(BaseField *p, int q){return p->velocityY[q];}
