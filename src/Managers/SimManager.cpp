#include <csignal>
#include <unistd.h>
#include "SimManager.h"
#include "../Backends/BackendFactory.h"
#include "../Utilities/RCexception.h"
#include "../Utilities/Timings.h"

void gbl_terminate(int arg) {
	// if the simulation has not started yet, then we make it so pressing ctrl+c twice
	// kills the program no matter what.
	if(!SimManager::started) signal(arg, SIG_DFL);
	OX_LOG(Logger::LOG_INFO, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	SimManager::stop = true;
}

bool SimManager::stop = false;
bool SimManager::started = false;

SimManager::SimManager(input_file inp) :
				input(inp),
				print_energy_every(1000) {
	steps = steps_run = equilibration_steps = 0;
	time_scale_manager.state = 0;
	print_input = 0;
	pid = getpid();
	seed = -1;
	backend = nullptr;
	time_scale_var = -1;
}

SimManager::~SimManager() {
	// print unread (i.e. unused) keys
	input.set_unread_keys();
	std::string unread;
	for(std::vector<std::string>::iterator it = input.unread_keys.begin(); it != input.unread_keys.end(); it++) {
		unread += std::string("\n\t") + *it;
	}
	if(unread.size() > 0) {
		OX_DEBUG("The following keys found in the input file were not used: %s", unread.c_str());
	}

	cleanTimeScale(&time_scale_manager);

	if(backend != nullptr) {
		int updated = backend->get_N_updates();
		if(updated > 0) {
			OX_LOG(Logger::LOG_INFO, "Lists updated %d times (every ~%lf steps)", updated, steps_run / (double)updated);
		}
		backend.reset();
	}
}

void SimManager::load_options() {
	Logger::instance()->get_settings(input);

	getInputBoolAsInt(&input, "print_input", &print_input, 0);
	getInputInt(&input, "print_energy_every", &print_energy_every, 0);
	getInputLLInt(&input, "steps", &steps, 1);
	getInputLLInt(&input, "equilibration_steps", &equilibration_steps, 0);

	// check that equilibration is only run on a simulation 
	bool my_restart_step_counter = false;
	getInputBool(&input, "restart_step_counter", &my_restart_step_counter, 0);
	if(my_restart_step_counter == false && equilibration_steps > 0) {
		throw RCexception("Incompatible key values found:\n\tif equilibration_steps > 0, restart_step_counter must be set to true.\n\tAborting");
	}

	if(equilibration_steps < 0) throw RCexception("Equilibration steps can not be < 0. Aborting");
	if(getInputInt(&input, "seed", &seed, 0) == KEY_NOT_FOUND) {
		seed = time(NULL);
		int rand_seed = 0;
		FILE *f = fopen("/dev/urandom", "rb");
		if(f == NULL) {
			OX_LOG(Logger::LOG_INFO, "Can't open /dev/urandom, using system time as a random seed");
		}
		else {
			if(fread((void *) &rand_seed, sizeof(rand_seed), 1, f) != 0) seed += rand_seed;
			else OX_LOG(Logger::LOG_INFO, "Can't read from /dev/urandom, using system time as a random seed");
			fclose(f);
		}
	}
}

void SimManager::init() {
	OX_LOG(Logger::LOG_INFO, "seeding the RNG with %d", seed);
	srand48(seed);

	OX_LOG(Logger::LOG_INFO, "Initializing backend ", seed);
	backend = BackendFactory::make_backend(input);
	backend->get_settings(input);

	// here we handle a few SIG* signals;
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	if(print_input) {
		char out_name[512];
		sprintf(out_name, "input.%d", pid);
		input.print(out_name);
	}
	OX_LOG(Logger::LOG_INFO, "pid: %d", pid);

	char ts_type[256];
	getInputString(&input, "time_scale", ts_type, 1);
	time_scale_var = TS_LIN;
	if(strcmp(ts_type, "linear") == 0) time_scale_var = TS_LIN;
	else if(strcmp(ts_type, "log_lin") == 0) time_scale_var = TS_LOG_LIN;
	else throw RCexception("Time scale '%s' not supported", ts_type);

	backend->init();

	// init time_scale_manager
	initTimeScale(&time_scale_manager, time_scale_var);

	int tmp, tmpm;
	getInputInt(&input, "print_conf_interval", &tmpm, 1);
	setTSInterval(&time_scale_manager, tmpm);

	if(time_scale_var == TS_LOG_LIN) {
		getInputInt(&input, "print_conf_ppc", &tmp, 1);
		setTSPPC(&time_scale_manager, tmp);
	}
	// the awkward second argument in the next line makes consistent
	// trajectory files...
	setTSInitialStep(&time_scale_manager, backend->start_step_from_file + tmpm - (backend->start_step_from_file % tmpm));
	// end
}

void SimManager::run() {
	backend->apply_changes_to_simulation_data();

	SimManager::started = true;
	// equilibration loop
	if(equilibration_steps > 0) {
		OX_LOG(Logger::LOG_INFO, "Equilibrating...");
		for(llint step = 0; step < equilibration_steps && !SimManager::stop; step++) {
			backend->sim_step();
		}
		OX_LOG(Logger::LOG_INFO, "Equilibration done");
	}
	backend->print_equilibration_info();
	backend->apply_simulation_changes_after_equilibration();

	// main loop
	OX_LOG(Logger::LOG_INFO, "Starting main loop...");
	for(steps_run = 0; steps_run < steps && !SimManager::stop; steps_run++) {
		if(backend->current_step() == time_scale_manager.next_step) {
		//if(backend->current_step()>5600 && backend->current_step()<5800) {
			if(steps_run > 0) {
				backend->print_conf();
			}
			setTSNextStep(&time_scale_manager);
		}

		backend->update_observables_data();
		backend->print_observables();
		backend->sim_step();
		backend->increment_current_step();
	}
	// this is in case _cur_step, after being increased by 1 before exiting the loop,
	// has become a multiple of print_conf_every
	if(backend->current_step() == time_scale_manager.next_step) {
		if(steps_run > 0) {
			backend->print_conf();
		}
		setTSNextStep(&time_scale_manager);
	}
	backend->print_observables();

	// prints the last configuration
	backend->print_conf(true);

	TimingManager::instance()->print(steps_run);

	backend->apply_simulation_data_changes();
}


void SimManager::analyse() {
	backend->apply_simulation_changes_after_equilibration();
	backend->update_observables_data();
	backend->print_observables();
}
