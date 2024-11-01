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
}

void FD_CPUBackend::init() {

	FDBackend::init();

	interaction->begin_energy_computation();
	//compute_forces();

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
	std::vector<int> particles_with_warning;

	for(auto p : fields) {
		//printf("first: %f, %f, %f, %f\n", p->area, p->sumF, p->CoM[0], p->CoM[1]);
		p->S00=0;
		p->S01=0;
		p->area=0;
		p->sumF=0;
		com_old[0]=p->CoM[0];
		com_old[1]=p->CoM[1];
		p->CoM[0] = 0; p->CoM[1] = 0;
		temp[0]=0; temp[1]=0; temp[2]=0; temp[3]=0;
		p->x_sub_left = box->getXsize();
		p->y_sub_bottom = box->getYsize();

		for(int q=0; q<p->subSize; q++) { 

			//int  k  = p->GetSubIndex(q, config_info->box);
			int k = p->map_sub_to_box[q];

			// compute dphi
	                dphi =
	                // free energy
	                - J0 * p->freeEnergy[q] 
               		// advection term
	                - interaction->get_velocity_x(p,q) * p->fieldDX[q] - interaction->get_velocity_y(p,q) * p->fieldDY[q];

			//if(p->index==2 && store==true)printf("check: %d, %d | %f, %f, %f, %f, %f | %f, %f, %d, %d\n", p->index, q, p->freeEnergy[q], p->velocityX, p->velocityY, p->fieldDX[q], p->fieldDY[q], com_old[0], com_old[1], p->offset[0],  p->offset[1]);

			if(store) {
				p->fieldScalar_old[q]=p->fieldScalar[q];
				p->dfield_old[q]=dphi;
			}

			phi = p->fieldScalar_old[q] + dt*.5*(dphi + p->dfield_old[q]);

			p->fieldScalar[q] = phi;
			p->area += phi * phi;
			p->sumF += phi;

			//timer_testing->resume();
			temp[0] += phi * cos_x_table[p->map_sub_to_box_x[q]];//GetSubXIndex(q, config_info->box)]; //cos(2*PI*p->GetSubXIndex(q, config_info->box)/box->getXsize()); 
			temp[1] += phi * sin_x_table[p->map_sub_to_box_x[q]];//GetSubXIndex(q, config_info->box)]; //sin(2*PI*p->GetSubXIndex(q, config_info->box)/box->getXsize());
			temp[2] += phi * cos_y_table[p->map_sub_to_box_y[q]];//GetSubYIndex(q, config_info->box)]; //cos(2*PI*p->GetSubYIndex(q, config_info->box)/box->getYsize());
			temp[3] += phi * sin_y_table[p->map_sub_to_box_y[q]];//GetSubYIndex(q, config_info->box)]; //sin(2*PI*p->GetSubYIndex(q, config_info->box)/box->getYsize());
			//timer_testing->pause();
			if(phi>0.1 && phi<0.2){
				dd = (p->sub_corner_bottom_left-int(p->sub_corner_bottom_left/box->getXsize())*box->getXsize()) - p->map_sub_to_box_x[q];
				if(dd>box->getXsize()/2)dd-=box->getXsize();
				else if(dd<-box->getXsize()/2)dd+=box->getXsize();
				if(abs(dd) < p->x_sub_left){p->x_sub_left=abs(dd);}

				dd = (p->sub_corner_bottom_left/box->getXsize()) - p->map_sub_to_box_y[q];
				if(dd>box->getYsize()/2)dd-=box->getYsize();
				else if(dd<-box->getYsize()/2)dd+=box->getYsize();
				if(abs(dd) < p->y_sub_bottom){p->y_sub_bottom=abs(dd);}
			}


			interaction->resetSums(k);
			interaction->updateFieldProperties(p, q, k);
		}
		p->CoM[0] = box->getXsize()*( atan2(-temp[1]/p->sumF, -temp[0]/p->sumF) + PI ) / (2*PI); 
		p->CoM[1] = box->getYsize()*( atan2(-temp[3]/p->sumF, -temp[2]/p->sumF) + PI ) / (2*PI);

		if(box->sqr_min_image_distance(com_old, p->CoM) > distance_thresh) {
			particles_with_warning.push_back(p->index);
		}

		//std::cout<< p->y_sub_bottom << " "<< p->x_sub_left << " "<< p->sub_corner_bottom_left/box->getXsize() << " "<<p->sub_corner_bottom_left-int(p->sub_corner_bottom_left/box->getXsize())*box->getXsize()<<" "<< p->map_sub_to_box_y[save_q]<<" "<<p->map_sub_to_box_x[save_q]<<std::endl;
		interaction->updateDirectedActiveForces(dt, p, store);
		if(!store)p->set_positions(config_info->box);
	}

	if(particles_with_warning.size() > 0) {
		std::stringstream ss;
		for(auto idx : particles_with_warning) {
			ss << idx << " ";
		}
		OX_LOG(Logger::LOG_WARNING, "The following particles had a displacement greater than %f in this step: %s", distance_thresh, ss.str().c_str());
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
