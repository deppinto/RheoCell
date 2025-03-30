#include <sstream>
#include <typeinfo>
#include <map>

#include "ForceField.h"
#include "../Fields/BaseField.h"
#include "../Fields/MultiPhaseField.h"
#include "../Fields/LEBcMultiPhaseField.h"

using namespace std;

ForceField::ForceField() {
	only_type = -1;
}

ForceField::~ForceField() {

}

void ForceField::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "only_type", &only_type, 0);
	getInputInt(&my_inp, "size_grid", &size_grid, 0);
}

void ForceField::init() {
	BaseObservable::init();

	Lx = config_info->box->getXsize();
	Ly = config_info->box->getYsize();
	f_field_x.resize(Lx * Ly);
	f_field_y.resize(Lx * Ly);
	f_field_passive_x.resize(Lx * Ly);
	f_field_passive_y.resize(Lx * Ly);
	f_field_active_x.resize(Lx * Ly);
	f_field_active_y.resize(Lx * Ly);
	phi_field.resize(Lx * Ly);


	int Lx_coarse = Lx; //int(Lx / delta)+1;
	int Ly_coarse = Ly; //int(Ly / delta)+1;
	f_field_coarse_x.resize(Lx_coarse * Ly_coarse);
	f_field_coarse_y.resize(Ly_coarse * Ly_coarse);
	f_field_passive_coarse_x.resize(Lx_coarse * Ly_coarse);
	f_field_passive_coarse_y.resize(Ly_coarse * Ly_coarse);
	f_field_active_coarse_x.resize(Lx_coarse * Ly_coarse);
	f_field_active_coarse_y.resize(Ly_coarse * Ly_coarse);
	phi_field_coarse.resize(Lx_coarse * Ly_coarse);

	if(visible_fields.size() == 0) {
                for(int i = 0; i < config_info->N(); i++) {
                        visible_fields.insert(i);
                }
        }
}

string ForceField::headers(llint step) {
	stringstream headers;
	headers.precision(15);

	headers << "t = " << step << endl;
	headers << "b = " << Lx << " " << Ly << endl;

	return headers.str();
}

string ForceField::field(BaseField *p) {
	stringstream conf;
	conf.precision(15);
	
	for(int i=0;i<p->subSize; i++)conf << p->map_sub_to_box[i] << " " << config_info->interaction->get_velocity_x(p, i) << " "<< config_info->interaction->get_velocity_y(p, i) <<" ";

	return conf.str();
}

void ForceField::calc_field(BaseField *p) {

	//std::cout<<"begin-------------------------------------------------------------------- "<<p->index<<std::endl;
	for(int q=0;q<p->subSize; q++){
		int k = p->map_sub_to_box[q];
		//std::cout<<"first: "<<k<<" "<<Lx<<" "<<Ly<<std::endl;

		//f_field_x[k] += (p->Fpassive_x[q] + p->Factive_x[q]) * p->fieldScalar[q];
		//f_field_y[k] += (p->Fpassive_y[q] + p->Factive_y[q]) * p->fieldScalar[q];
		f_field_x[k] += config_info->interaction->get_total_force_x(p, q) * p->fieldScalar[q];
		f_field_y[k] += config_info->interaction->get_total_force_y(p, q) * p->fieldScalar[q];
		f_field_passive_x[k] += p->Fpassive_x[q] * p->fieldScalar[q];
		f_field_passive_y[k] += p->Fpassive_y[q] * p->fieldScalar[q];
		f_field_active_x[k] += p->Factive_x[q] * p->fieldScalar[q];
		f_field_active_y[k] += p->Factive_y[q] * p->fieldScalar[q];
		phi_field[k] += p->fieldScalar[q];

		//int y = int(p->map_sub_to_box_y[q]/delta);
		//int x = int(p->map_sub_to_box_x[q]/delta);
		//k = x + y * Lx_coarse;
		//std::cout<<"second: "<<k<<" "<<Lx_coarse<<" "<<Ly_coarse<<std::endl;

		/*f_field_coarse_x[k] += (p->Fpassive_x[q] + p->Factive_x[q]) * p->fieldScalar[q];
		f_field_coarse_y[k] += (p->Fpassive_y[q] + p->Factive_y[q]) * p->fieldScalar[q];
		f_field_passive_coarse_x[k] += p->Fpassive_x[q] * p->fieldScalar[q];
		f_field_passive_coarse_y[k] += p->Fpassive_y[q] * p->fieldScalar[q];
		f_field_active_coarse_x[k] += p->Factive_x[q] * p->fieldScalar[q];
		f_field_active_coarse_y[k] += p->Factive_y[q] * p->fieldScalar[q];
		phi_field_coarse[k] += p->fieldScalar[q];*/
	}
	//std::cout<<"cell-------------------------------------------------------------------- "<<p->index<<std::endl;
}

string ForceField::f_field(llint step) {
	stringstream conf;
	conf.precision(15);
	std::fill(f_field_x.begin(), f_field_x.end(), 0.);
	std::fill(f_field_y.begin(), f_field_y.end(), 0.);
	std::fill(phi_field.begin(), phi_field.end(), 0.);
	std::fill(f_field_coarse_x.begin(), f_field_coarse_x.end(), 0.);
	std::fill(f_field_coarse_y.begin(), f_field_coarse_y.end(), 0.);
	std::fill(phi_field_coarse.begin(), phi_field_coarse.end(), 0.);
	std::fill(f_field_passive_x.begin(), f_field_passive_x.end(), 0.);
	std::fill(f_field_passive_y.begin(), f_field_passive_y.end(), 0.);
	std::fill(f_field_active_x.begin(), f_field_active_x.end(), 0.);
	std::fill(f_field_active_y.begin(), f_field_active_y.end(), 0.);
	std::fill(f_field_passive_coarse_x.begin(), f_field_passive_coarse_x.end(), 0.);
	std::fill(f_field_passive_coarse_y.begin(), f_field_passive_coarse_y.end(), 0.);
	std::fill(f_field_active_coarse_x.begin(), f_field_active_coarse_x.end(), 0.);
	std::fill(f_field_active_coarse_y.begin(), f_field_active_coarse_y.end(), 0.);

	for(auto p_idx : visible_fields) {
		BaseField *p = config_info->fields()[p_idx];
		bool visible = (only_type == -1 || p->type == only_type);
		if(visible) {
			calc_field(p);
		}
	}

	for(int i=0; i<Ly; i++){
		for(int j=0; j<Lx; j++){
			int site = j + i * Lx;
			conf << j+0.5 <<" "<< i+0.5 << " " << f_field_x[site] << " " << f_field_y[site] << " "<< " "<< f_field_passive_x[site]<<" "<<f_field_passive_y[site]<<" "<<f_field_active_x[site]<<" "<<f_field_active_y[site]<<" ";
		}
	}

	conf<<"\n";
	for(int i=0; i<Ly; i++){
		for(int j=0; j<Lx; j++){
			int site = j + i * Lx;

			for(int k=0; k<size_grid; k++){
				for(int l=0; l<size_grid; l++){

					int xx = (((j - int(size_grid/2) + Lx) % Lx) + l) % Lx;
					int yy = (((i - int(size_grid/2) + Ly) % Ly) + k) % Ly;
					int ss = xx + yy * Lx;

					f_field_coarse_x[site] += (f_field_x[ss]) / (size_grid * size_grid);
					f_field_coarse_y[site] += (f_field_y[ss]) /(size_grid * size_grid) ;
					f_field_passive_coarse_x[site] += (f_field_passive_x[ss]) /(size_grid * size_grid) ;
					f_field_passive_coarse_y[site] += (f_field_passive_y[ss]) /(size_grid * size_grid) ;
					f_field_active_coarse_x[site] += (f_field_active_x[ss]) /(size_grid * size_grid) ;
					f_field_active_coarse_y[site] += (f_field_active_y[ss]) /(size_grid * size_grid) ;
					phi_field_coarse[site] += 1.;				
				}
			}

			conf << j + 0.5 <<" "<<  i + 0.5  << " " << f_field_coarse_x[site] << " " << f_field_coarse_y[site] << " "<< f_field_passive_coarse_x[site]<<" "<<f_field_passive_coarse_y[site]<<" "<<f_field_active_coarse_x[site]<<" "<<f_field_active_coarse_y[site]<<" ";
		}
	}

	return conf.str();
}

string ForceField::get_output_string(llint curr_step) {
	return headers(curr_step) + f_field(curr_step);
}
