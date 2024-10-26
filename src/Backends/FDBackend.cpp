#include <sstream>

#include "FDBackend.h"
#include "../Observables/ObservableOutput.h"
#include "../Utilities/Timings.h"


FDBackend::FDBackend() :
				SimBackend() {
	timer_first_step = timer_forces = timer_testing = NULL;
	lees_edwards = false;
	shear_rate = -0.f;
}

FDBackend::~FDBackend() {

}

void FDBackend::get_settings(input_file &inp) {
	SimBackend::get_settings(inp);

	getInputBool(&inp, "lees_edwards", &lees_edwards, 0);
	if(lees_edwards) {
		getInputNumber(&inp, "lees_edwards_shear_rate", &shear_rate, 1);
		if(shear_rate <= 0.) throw RCexception("lees_edwards_shear_rate should be > 0");
		OX_LOG(Logger::LOG_INFO, "Using Lees-Edwards boundary conditions with shear rate %lf", shear_rate);
	}

	getInputNumber(&inp, "dt", &dt, 1);
}


void FDBackend::init() {
	SimBackend::init();

	timer_first_step = TimingManager::instance()->new_timer(string("First Step"), string("SimBackend"));
	timer_forces = TimingManager::instance()->new_timer(string("Forces"), string("SimBackend"));
	timer_testing = TimingManager::instance()->new_timer(string("Testing"), string("SimBackend"));
}

void FDBackend::print_observables() {
	SimBackend::print_observables();
}
