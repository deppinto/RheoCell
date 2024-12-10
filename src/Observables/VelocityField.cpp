#include <sstream>
#include <typeinfo>
#include <map>

#include "VelocityField.h"
#include "../Fields/BaseField.h"
#include "../Fields/MultiPhaseField.h"
#include "../Fields/LEBcMultiPhaseField.h"

using namespace std;

VelocityField::VelocityField() {
	only_type = -1;
}

VelocityField::~VelocityField() {

}

void VelocityField::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "only_type", &only_type, 0);
}

void VelocityField::init() {
	BaseObservable::init();

	Lx = config_info->box->getXsize();
	Ly = config_info->box->getYsize();
	v_field_x.resize(Lx * Ly);
	v_field_y.resize(Lx * Ly);
	phi_field.resize(Lx * Ly);

	if(visible_fields.size() == 0) {
                for(int i = 0; i < config_info->N(); i++) {
                        visible_fields.insert(i);
                }
        }
}

string VelocityField::headers(llint step) {
	stringstream headers;
	headers.precision(15);

	headers << "t = " << step << endl;
	headers << "b = " << Lx << " " << Ly << endl;

	return headers.str();
}

string VelocityField::field(BaseField *p) {
	stringstream conf;
	conf.precision(15);
	
	for(int i=0;i<p->subSize; i++)conf << p->GetSubIndex(i, config_info->box) << " " << config_info->interaction->get_velocity_x(p, i) << " "<< config_info->interaction->get_velocity_y(p, i) <<" ";

	return conf.str();
}

void VelocityField::calc_velocity_field(BaseField *p) {

	for(int i=0;i<p->subSize; i++){
		number velocityX = config_info->interaction->get_velocity_x(p, i);
		number velocityY = config_info->interaction->get_velocity_y(p, i);
		v_field_x[p->GetSubIndex(i, config_info->box)] += velocityX; //* p->fieldScalar[i];
		v_field_y[p->GetSubIndex(i, config_info->box)] += velocityY; //* p->fieldScalar[i];
		phi_field[p->GetSubIndex(i, config_info->box)] += p->fieldScalar[i];
	}
}

string VelocityField::velocity_field(llint step) {
	stringstream conf;
	conf.precision(15);
	std::fill(v_field_x.begin(), v_field_x.end(), 0.);
	std::fill(v_field_y.begin(), v_field_y.end(), 0.);
	std::fill(phi_field.begin(), phi_field.end(), 0.);

	for(auto p_idx : visible_fields) {
		BaseField *p = config_info->fields()[p_idx];
		bool visible = (only_type == -1 || p->type == only_type);
		if(visible) {
			calc_velocity_field(p);
		}
	}

	for(int i=0; i<Ly; i++){
		for(int j=0; j<Lx; j++){
			int site = j + i * Lx;
			if(phi_field[site]>0)
				conf << site << " " << v_field_x[site] << " " << v_field_y[site] << " ";
				//conf << site << " " << v_field_x[site]/phi_field[site] << " " << v_field_y[site]/phi_field[site] << " ";
			else
				conf << site << " " << 0 << " " << 0 << " ";
		}
	}

	return conf.str();
}

string VelocityField::get_output_string(llint curr_step) {
	return headers(curr_step) + velocity_field(curr_step);
}