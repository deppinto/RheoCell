#include <sstream>
#include "FD_CPUBackend.h"
#include "../Utilities/Timings.h"

FD_CPUBackend::FD_CPUBackend() :
				FDBackend() {
}

FD_CPUBackend::~FD_CPUBackend() {

}

void FD_CPUBackend::get_settings(input_file &inp) {
	FDBackend::get_settings(inp);
	getInputNumber(&inp, "J0", &J0, 0);
	getInputBool(&inp, "analysis", &analysis, 0);
}

void FD_CPUBackend::init() {

	FDBackend::init();

	interaction->begin_energy_computation();
	if(analysis){compute_forces();}

	int sizeX = box->getXsize();
	int sizeY = box->getYsize();
	cos_x_table.resize(sizeX);
	cos_y_table.resize(sizeY);
	sin_x_table.resize(sizeX);
	sin_y_table.resize(sizeY);
	for(int i =0; i<sizeX; i++){
		cos_x_table[i]=cos(2*PI*i/sizeX);
		sin_x_table[i]=sin(2*PI*i/sizeX);
	}
        for(int i =0; i<sizeY; i++){
                cos_y_table[i]=cos(2*PI*i/sizeY);
                sin_y_table[i]=sin(2*PI*i/sizeY);
        }

}

void FD_CPUBackend::first_step(bool store) {

	for(auto p : fields) {
		com_old[0]=p->CoM[0];
		com_old[1]=p->CoM[1];
		temp[0]=0; temp[1]=0; temp[2]=0; temp[3]=0;
		number area = p->area;
		number sumF = p->sumF;
		p->set_properties_to_zero();

		//number sumTest_x=0.;
		//number sum_x=0.;

		for(int q=0; q<p->subSize; q++) { 

			//int  k  = p->GetSubIndex(q, config_info->box);
			int k = p->map_sub_to_box[q];

			// compute dphi
	                dphi =
	                // free energy
	                - J0 * p->freeEnergy[q] 
               		// advection term
	                - interaction->get_velocity_x(p,q) * p->fieldDX[q] - interaction->get_velocity_y(p,q) * p->fieldDY[q];

			if(store) {
				p->fieldScalar_old[q]=p->fieldScalar[q];
				p->dfield_old[q]=dphi;
				p->area_old = area;
				p->sumF_old = sumF;
			}

			phi = p->fieldScalar_old[q] + dt*.5*(dphi + p->dfield_old[q]);

			//if(int(q/p->LsubX)>3){sumTest_x+=phi * p->map_sub_to_box_x[q];sum_x+=phi;}

			p->fieldScalar[q] = phi;
			p->area += phi * phi;
			p->sumF += phi;
			
			//timer_testing->resume();
			int x=(((q - int(q/p->LsubX) * p->LsubX)+p->offset[0])%p->LsubX);
			int y=(((q/p->LsubX)+p->offset[1])%p->LsubY);
			temp[0] += phi * p->cos_x_table[x];//cos_x_table[p->map_sub_to_box_x[q]];//GetSubXIndex(q, config_info->box)]; //cos(2*PI*p->GetSubXIndex(q, config_info->box)/box->getXsize()); 
			temp[1] += phi * p->sin_x_table[x];//sin_x_table[p->map_sub_to_box_x[q]];//GetSubXIndex(q, config_info->box)]; //sin(2*PI*p->GetSubXIndex(q, config_info->box)/box->getXsize());
			temp[2] += phi * p->cos_y_table[y];//cos_y_table[p->map_sub_to_box_y[q]];//GetSubYIndex(q, config_info->box)]; //cos(2*PI*p->GetSubYIndex(q, config_info->box)/box->getYsize());
			temp[3] += phi * p->sin_y_table[y];//sin_y_table[p->map_sub_to_box_y[q]];//GetSubYIndex(q, config_info->box)]; //sin(2*PI*p->GetSubYIndex(q, config_info->box)/box->getYsize());
			//timer_testing->pause();

			if(!store && phi<0.2 && phi>0.1)p->check_borders(q);
			interaction->resetSums(k);
			interaction->updateFieldProperties(p, q, k);
		}
		//p->CoM[0] = box->getXsize()*( atan2(-temp[1]/p->sumF, -temp[0]/p->sumF) + PI ) / (2*PI); 
		//p->CoM[1] = box->getYsize()*( atan2(-temp[3]/p->sumF, -temp[2]/p->sumF) + PI ) / (2*PI);
		p->CoM[0] = p->LsubX*( atan2(-temp[1]/p->sumF, -temp[0]/p->sumF) + PI ) / (2*PI); 
		p->CoM[1] = p->LsubY*( atan2(-temp[3]/p->sumF, -temp[2]/p->sumF) + PI ) / (2*PI);


		//std::cout<<"backend: "<<p->fieldScalar[p->subSize-1]<<" "<<p->freeEnergy[p->subSize-1]<<" "<<interaction->get_velocity_x(p,p->subSize-1) <<" "<<interaction->get_velocity_y(p,p->subSize-1) <<" "<<p->fieldDX[p->subSize-1]<<" "<<p->fieldDY[p->subSize-1]<<" "<<p->CoM[0]<<" "<<p->CoM[1]<<" "<<config_info->curr_step  <<std::endl;
		std::cout<<"backend: "<<p->fieldScalar[0]<<" "<<p->freeEnergy[0]<<" "<<interaction->get_velocity_x(p,0) <<" "<<interaction->get_velocity_y(p,0) <<" "<<p->fieldDX[0]<<" "<<p->fieldDY[0]<<" "<<p->CoM[0]<<" "<<p->CoM[1]<<" "<<config_info->curr_step << " " << p->index << " "<<p->sub_corner_bottom_left/box->getXsize() <<std::endl;

		interaction->updateDirectedActiveForces(dt, p, store);
		if(!store)p->set_positions(config_info->box);
		std::cout<<"done"<<std::endl;
	}
}

void FD_CPUBackend::compute_forces() {

	interaction->begin_energy_computation(fields);
	U = interaction->U;
	K = interaction->K;
}


void FD_CPUBackend::second_step() {

}

void FD_CPUBackend::sim_step() {
	mytimer->resume();

	timer_forces->resume();
	compute_forces();
	second_step();
	timer_forces->pause();

	timer_first_step->resume();
	first_step(true);
	timer_first_step->pause();

        timer_forces->resume();
        compute_forces();
        second_step();
        timer_forces->pause();

        timer_first_step->resume();
        first_step(false);
        timer_first_step->pause();

	mytimer->pause();
}
