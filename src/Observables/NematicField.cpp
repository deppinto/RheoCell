#include <sstream>
#include <typeinfo>
#include <map>

#include "NematicField.h"
#include "../Fields/BaseField.h"
#include "../Fields/MultiPhaseField.h"
#include "../Fields/LEBcMultiPhaseField.h"

using namespace std;

NematicField::NematicField() {
	only_type = -1;
}

NematicField::~NematicField() {

}

void NematicField::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "only_type", &only_type, 0);
	getInputInt(&my_inp, "size_grid", &size_grid, 0);
}

void NematicField::init() {
	BaseObservable::init();

	Lx = config_info->box->getXsize();
	Ly = config_info->box->getYsize();
	f_field_x.resize(Lx * Ly);
	f_field_y.resize(Lx * Ly);
	phi_field.resize(Lx * Ly);


	int Lx_coarse = Lx; //int(Lx / delta)+1;
	int Ly_coarse = Ly; //int(Ly / delta)+1;
	f_field_coarse_x.resize(Lx_coarse * Ly_coarse);
	f_field_coarse_y.resize(Ly_coarse * Ly_coarse);
	phi_field_coarse.resize(Lx_coarse * Ly_coarse);

	if(visible_fields.size() == 0) {
                for(int i = 0; i < config_info->N(); i++) {
                        visible_fields.insert(i);
                }
        }
}

string NematicField::headers(llint step) {
	stringstream headers;
	headers.precision(15);

	headers << "t = " << step << endl;
	headers << "b = " << Lx << " " << Ly << endl;

	return headers.str();
}

string NematicField::field(BaseField *p) {
	stringstream conf;
	conf.precision(15);
	
	for(int i=0;i<p->subSize; i++)conf << p->map_sub_to_box[i] << " " << config_info->interaction->get_velocity_x(p, i) << " "<< config_info->interaction->get_velocity_y(p, i) <<" ";

	return conf.str();
}

void NematicField::calc_field(BaseField *p) {

	for(int q=0;q<p->subSize; q++){
		int k = p->map_sub_to_box[q];

		f_field_x[k] += p->Q00 * p->fieldScalar[q];
		f_field_y[k] += p->Q01 * p->fieldScalar[q];
		phi_field[k] += p->fieldScalar[q];
	}
}

string NematicField::f_field(llint step) {
	stringstream conf;
	conf.precision(15);
	std::fill(f_field_x.begin(), f_field_x.end(), 0.);
	std::fill(f_field_y.begin(), f_field_y.end(), 0.);
	std::fill(phi_field.begin(), phi_field.end(), 0.);
	std::fill(f_field_coarse_x.begin(), f_field_coarse_x.end(), 0.);
	std::fill(f_field_coarse_y.begin(), f_field_coarse_y.end(), 0.);
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
			conf << j+0.5 <<" "<< i+0.5 << " " << f_field_x[site] << " " << f_field_y[site] << " ";
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

					f_field_coarse_x[site] += f_field_x[ss] / (size_grid * size_grid);
					f_field_coarse_y[site] += f_field_y[ss] / (size_grid * size_grid);
					phi_field_coarse[site] += 1.;				
				}
			}

			conf << j + 0.5 <<" "<<  i + 0.5  << " " << f_field_coarse_x[site] << " " << f_field_coarse_y[site] << " ";
		}
	}

	return conf.str();
}

string NematicField::get_output_string(llint curr_step) {
	return headers(curr_step) + f_field(curr_step);
}
