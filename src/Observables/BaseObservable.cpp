/*
 * BaseObservable.cpp
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#include "BaseObservable.h"

BaseObservable::BaseObservable() :
				config_info(ConfigInfo::instance().get()) {
}

BaseObservable::~BaseObservable() {

}

bool BaseObservable::need_updating(llint curr_step) {
	return (int_update_every > 0 && (curr_step % int_update_every) == 0);
}

void BaseObservable::update_data(llint curr_step) {

}

void BaseObservable::get_settings(input_file &my_inp, input_file &sim_inp) {
	getInputString(&my_inp, "id", id, 0);
	getInputLLInt(&my_inp, "update_every", &int_update_every, 0);
}

void BaseObservable::init() {

}
