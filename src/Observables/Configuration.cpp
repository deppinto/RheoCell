#include <sstream>
#include <typeinfo>
#include <map>

#include "Configuration.h"
#include "../Fields/BaseField.h"
#include "../Fields/MultiPhaseField.h"

using namespace std;

Configuration::Configuration() {
	only_type = -1;
}

Configuration::~Configuration() {

}

void Configuration::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "only_type", &only_type, 0);
}

void Configuration::init() {
	BaseObservable::init();

	if(visible_fields.size() == 0) {
                for(int i = 0; i < config_info->N(); i++) {
                        visible_fields.insert(i);
                }
        }
}

string Configuration::headers(llint step) {
	stringstream headers;
	headers.precision(15);

	headers << "t = " << step << endl;
	headers << "b = " << config_info->box->getXsize() << " " << config_info->box->getYsize() << endl;

	return headers.str();
}

string Configuration::field(BaseField *p) {
	stringstream conf;
	conf.precision(15);
	
	conf << p->LsubX<< " "<< p->LsubY<< " "<<p->CoM[0] << " " << p->CoM[1] << " " << p->offset[0] << " " << p->offset[1] << " " << p->sub_corner_bottom_left << " " << p->nemQ[0] << " "<< p->nemQ[1] << " ";
	for(int i=0;i<p->subSize; i++)conf << p->GetSubIndex(i, config_info->box) << " " << p->fieldScalar[i] << " ";
	//conf << endl;

	return conf.str();
}

string Configuration::configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	// this is used to avoid printing empty lines
	bool empty = true;
	for(auto p_idx : visible_fields) {
		BaseField *p = config_info->fields()[p_idx];
		bool visible = (only_type == -1 || p->type == only_type);
		if(visible) {
			if(p_idx != *visible_fields.begin() && !empty) conf << endl;
			string p_str = field(p);
			conf << p_str;
			empty = (p_str.size() == 0);
		}
	}

	return conf.str();
}

string Configuration::get_output_string(llint curr_step) {
	return headers(curr_step) + configuration(curr_step);
}
