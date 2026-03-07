#include <sstream>
#include <typeinfo>
#include <map>

#include "ThermodynamicStressField.h"
#include "../Fields/BaseField.h"
#include "../Fields/MultiPhaseField.h"
#include "../Fields/LEBcMultiPhaseField.h"

using namespace std;

ThermodynamicStressField::ThermodynamicStressField() {
	only_type = -1;
}

ThermodynamicStressField::~ThermodynamicStressField() {

}

void ThermodynamicStressField::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "only_type", &only_type, 0);
	getInputInt(&my_inp, "size_grid", &size_grid, 0);
}

void ThermodynamicStressField::init() {
	BaseObservable::init();

	Lx = config_info->box->getXsize();
	Ly = config_info->box->getYsize();
	f_field_xx.resize(Lx * Ly);
	f_field_yy.resize(Lx * Ly);
	f_field_xy.resize(Lx * Ly);
	phi_field.resize(Lx * Ly);


	int Lx_coarse = Lx; //int(Lx / delta)+1;
	int Ly_coarse = Ly; //int(Ly / delta)+1;
	f_field_coarse_xx.resize(Lx_coarse * Ly_coarse);
	f_field_coarse_yy.resize(Ly_coarse * Ly_coarse);
	f_field_coarse_xy.resize(Ly_coarse * Ly_coarse);
	phi_field_coarse.resize(Lx_coarse * Ly_coarse);

	if(visible_fields.size() == 0) {
                for(int i = 0; i < config_info->N(); i++) {
                        visible_fields.insert(i);
                }
        }
}

string ThermodynamicStressField::headers(llint step) {
	stringstream headers;
	headers.precision(15);

	headers << "t = " << step << endl;
	headers << "b = " << Lx << " " << Ly << endl;

	return headers.str();
}

string ThermodynamicStressField::field(BaseField *p) {
	stringstream conf;
	conf.precision(15);
	
	for(int i=0;i<p->subSize; i++)conf << p->map_sub_to_box[i] << " " << config_info->interaction->get_velocity_x(p, i) << " "<< config_info->interaction->get_velocity_y(p, i) <<" ";

	return conf.str();
}

void ThermodynamicStressField::calc_field(BaseField *p) {

	//std::cout<<"begin-------------------------------------------------------------------- "<<p->index<<std::endl;
	//number value_calc = 0.;
	//number value_calc_sides = 0.;
	//number value_calc_ends = 0.;
	for(int q=0;q<p->subSize; q++){
		//int xx = q%30;
		//int yy = q/30;
		int k = p->map_sub_to_box[q];
		//std::cout<<"first: "<<k<<" "<<Lx<<" "<<Ly<<std::endl;

		f_field_xx[k] += p->freeEnergyDensity[q] - p->freeEnergy[q] * p->fieldScalar[q] + p->fieldDX[q] * p->freeEnergyDensityGradient_x[q];
		f_field_yy[k] += p->freeEnergyDensity[q] - p->freeEnergy[q] * p->fieldScalar[q] + p->fieldDY[q] * p->freeEnergyDensityGradient_y[q];
		f_field_xy[k] += p->fieldDX[q] * p->freeEnergyDensityGradient_y[q];
		//value_calc += p->fieldDX[q] * p->freeEnergyDensityGradient_y[q];
		//if(xx<15 && yy>15)value_calc_sides += p->fieldDX[q] * p->freeEnergyDensityGradient_y[q];
		//else if(xx>15 && yy<15)value_calc_sides += p->fieldDX[q] * p->freeEnergyDensityGradient_y[q];
		//else value_calc_ends += p->fieldDX[q] * p->freeEnergyDensityGradient_y[q];
		phi_field[k] += p->fieldScalar[q];
	}
	//std::cout<<"cell-------------------------------------------------------------------- "<<p->index<<" "<<value_calc<<" "<<value_calc_sides<<" "<<value_calc_ends<<std::endl;
}

string ThermodynamicStressField::f_field(llint step) {
	stringstream conf;
	conf.precision(15);
	std::fill(f_field_xx.begin(), f_field_xx.end(), 0.);
	std::fill(f_field_yy.begin(), f_field_yy.end(), 0.);
	std::fill(f_field_xy.begin(), f_field_xy.end(), 0.);
	std::fill(phi_field.begin(), phi_field.end(), 0.);
	std::fill(f_field_coarse_xx.begin(), f_field_coarse_xx.end(), 0.);
	std::fill(f_field_coarse_yy.begin(), f_field_coarse_yy.end(), 0.);
	std::fill(f_field_coarse_xy.begin(), f_field_coarse_xy.end(), 0.);
	std::fill(phi_field_coarse.begin(), phi_field_coarse.end(), 0.);


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
			conf << j+0.5 <<" "<< i+0.5 << " " << f_field_xx[site] << " " << f_field_yy[site] << " "<< f_field_xy[site] <<" ";
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

					f_field_coarse_xx[site] += f_field_xx[ss] / (size_grid * size_grid);
					f_field_coarse_yy[site] += f_field_yy[ss] / (size_grid * size_grid);
					f_field_coarse_xy[site] += f_field_xy[ss] / (size_grid * size_grid);
					phi_field_coarse[site] += 1.;				
				}
			}


			conf << j + 0.5 <<" "<<  i + 0.5  << " " << f_field_coarse_xx[site] << " " << f_field_coarse_yy[site] << " "<< f_field_coarse_xy[site] << " ";
		}
	}

	return conf.str();
}


string ThermodynamicStressField::get_output_string(llint curr_step) {
	return headers(curr_step) + f_field(curr_step);
}
