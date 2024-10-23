#include <sstream>
#include <fstream>

#include "SimBackend.h"
#include "../Utilities/Utils.h"
#include "../Utilities/ConfigInfo.h"
#include "../Interactions/InteractionFactory.h"
#include "../Observables/ObservableFactory.h"
#include "../Forces/ForceFactory.h"
#include "../Boxes/BoxFactory.h"
#include "../Fields/BaseField.h"
#include "../Utilities/Timings.h"

SimBackend::SimBackend() {
	// we need to initialize everything so that we can check what we can
	// and what we can't delete[] in the destructor
	N_updates = 0;
	confs_to_skip = 0;
	interaction = nullptr;
	start_step_from_file = (llint) 0;
	backend_info = std::string("");
	reseed = false;
	custom_conf_name = false;
	read_conf_step = 0;
	box = nullptr;
	max_io = 1.e30;
	read_conf_step = -1;
	conf_interval = -1;
	obs_output_trajectory = obs_output_stdout = obs_output_file = obs_output_last_conf = nullptr;
	mytimer = nullptr;
	restart_step_counter = false;

	ConfigInfo::init(&fields);
	config_info = ConfigInfo::instance().get();
}

SimBackend::~SimBackend() {
	for(auto p : fields) {
		delete p;
	}

	// here we print the input output information
	llint total_file = 0;
	llint total_stderr = 0;
	OX_LOG(Logger::LOG_INFO, "Aggregated I/O statistics (set debug=1 for file-wise information)");
	for(auto const &element : obs_outputs) {
		llint now = element->get_bytes_written();
		auto fname = element->get_output_name();
		if(!strcmp(fname.c_str(), "stderr") || !strcmp(fname.c_str(), "stdout")) {
			total_stderr += now;
		}
		else {
			total_file += now;
		}
		std::string mybytes = Utils::bytes_to_human(now);
		OX_DEBUG("  on %s: %s", fname.c_str(), mybytes.c_str());
	}
	OX_LOG(Logger::LOG_NOTHING, "\t%s written to files" , Utils::bytes_to_human(total_file).c_str());
	OX_LOG(Logger::LOG_NOTHING, "\t%s written to stdout/stderr" , Utils::bytes_to_human(total_stderr).c_str());

	if(mytimer != nullptr) {
		double time_passed = (double) mytimer->get_time() / (double) CLOCKS_PER_SEC;
		if(time_passed > 0.) {
			OX_LOG(Logger::LOG_NOTHING, "\tFor a total of %8.3lg MB/s\n", (total_file + total_stderr) / ((1024.*1024.) * time_passed));
		}
	}
}

void SimBackend::get_settings(input_file &inp) {

	config_info->sim_input = &inp;

	// initialise the timer
	mytimer = TimingManager::instance()->new_timer(std::string("SimBackend"));

	interaction = InteractionFactory::make_interaction(inp);
	interaction->get_settings(inp);

	box = BoxFactory::make_box(inp);
	box->get_settings(inp);

	getInputBool(&inp, "restart_step_counter", &restart_step_counter, 0);

	// reload configuration
	std::string reload_from;
	if(getInputString(&inp, "reload_from", reload_from, 0) == KEY_FOUND) {

		// check that conf_file is not specified
		std::string tmpstring;
		if(getInputString(&inp, "conf_file", tmpstring, 0) == KEY_FOUND)
			throw RCexception("Input file error: \"conf_file\" cannot be specified if \"reload_from\" is specified");

		// check that restart_step_counter is set to 0
		if(restart_step_counter)
			throw RCexception("Input file error: \"restart_step_counter\" must be set to false if \"reload_from\" is specified");

		int my_seed;
		if(getInputInt(&inp, "seed", &my_seed, 0) == KEY_FOUND)
			throw RCexception("Input file error: \"seed\" must not be specified if \"reload_from\" is specified");

		conf_filename = std::string(reload_from);
	}
	else {
		getInputString(&inp, "conf_file", conf_filename, 1);
	}

	// we only reseed the RNG if:
	// a) we have a binary conf
	// b) no seed was specified in the input file
	// c) restart_step_counter is set to one (true) -- LR: looking at the code, it looks like it is quite the contrary
	int tmpi;
	getInputBoolAsInt(&inp, "restart_step_counter", &tmpi, 1);
	reseed = (getInputInt(&inp, "seed", &tmpi, 0) == KEY_NOT_FOUND || tmpi == 0);

	getInputInt(&inp, "confs_to_skip", &confs_to_skip, 0);

	char raw_T[256];
	getInputString(&inp, "T", raw_T, 1);
	config_info->update_temperature(Utils::get_temperature(raw_T));

	// here we fill the _obs_outputs vector
	auto new_outputs = ObservableFactory::make_observables();
	for(auto new_output : new_outputs) {
		add_output(new_output);
	}

	// we build the default stream of observables for trajectory and last configuration
	bool traj_print_momenta = false;
	getInputBool(&inp, "trajectory_print_momenta", &traj_print_momenta, 0);
	std::string traj_file;
	// Trajectory
	getInputString(&inp, "trajectory_file", traj_file, 1);
	std::string output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n}\n", traj_file.c_str());
	obs_output_trajectory = std::make_shared<ObservableOutput>(output_inp_text);

	std::string obs_text = Utils::sformat("type = configuration\nprint_momenta = %d", traj_print_momenta);
	obs_output_trajectory->add_observable(obs_text);
        add_output(obs_output_trajectory);

	// Last configuration
	std::string lastconf_file = "last_conf.dat";
	getInputString(&inp, "lastconf_file", lastconf_file, 0);
	output_inp_text = Utils::sformat("{\n\tname = %s\n\tprint_every = 0\n\tonly_last = 1\n}\n", lastconf_file.c_str());
	obs_output_last_conf = std::make_shared<ObservableOutput>(output_inp_text);
	obs_output_last_conf->add_observable("type = configuration\nid = last_conf");
	add_output(obs_output_last_conf);

	// set the max IO
	if(getInputNumber(&inp, "max_io", &max_io, 0) == KEY_FOUND) {
		if(max_io < 0)
			throw RCexception("Cannot run with a negative I/O limit. Set the max_io key to something > 0");
		else
			OX_LOG(Logger::LOG_INFO, "Setting the maximum IO limit to %g MB/s", max_io);
	}
	else {
		max_io = 1.; // default value for a simulation is 1 MB/s;
	}

}

void SimBackend::init() {
	conf_input.open(conf_filename.c_str());
	if(conf_input.good() == false) {
		throw RCexception("Can't read configuration file '%s'", conf_filename.c_str());
	}

	interaction->init();

	// check number of particles
	int N = interaction->get_N_from_topology();
	fields.resize(N);
	interaction->read_topology(fields);

	conf_input.seekg(0, std::ios::beg);

	// we need to skip a certain number of lines, depending on how many
	// particles we have and how many configurations we want to skip
	if(confs_to_skip > 0) {
		OX_LOG(Logger::LOG_INFO, "Skipping %d configuration(s)", confs_to_skip);
		int i;
		for(i = 0; i < confs_to_skip && conf_input.good(); i++) {
			if(!read_next_configuration(false)) {
				throw RCexception("Skipping %d configuration(s) is not possible, as the initial trajectory file only contains %d configurations", confs_to_skip, i);
			}
		}
	}

	bool check = read_next_configuration(false);
        if(!check) {
                throw RCexception("Could not read the initial configuration, aborting");
        }

	start_step_from_file = (restart_step_counter) ? 0 : read_conf_step;
	config_info->curr_step = start_step_from_file;

	// initialise external forces
	ForceFactory::instance()->make_forces(fields, box.get());

	interaction->set_box(box.get());

	config_info->set(interaction.get(), &backend_info, box.get());

	// initializes the observable output machinery. This part has to follow read_topology() since _fields has to be initialized
	for(auto const &element : obs_outputs) {
		element->init();
	}

}

bool SimBackend::read_next_configuration(bool binary) {
	double Lx, Ly;
	// parse headers. Binary and ascii configurations have different headers, and hence
	// we have to separate the two procedures
	{
		bool malformed_headers = false;
		std::string line;
		std::getline(conf_input, line);
		if(conf_input.eof()) {
			return false;
		}
		int res = sscanf(line.c_str(), "t = %lld", &read_conf_step);
		std::stringstream error_message;
		// handle the case when t can't be read
		if(res != 1) {
			error_message << "Malformed headers found in an input configuration.\"t = <int>\" was expected, but \"" << line << "\" was found instead.";
			malformed_headers = true;

		}
		if(!malformed_headers) {
			std::getline(conf_input, line);
			// handle the case when the box_size can't be read
			res += sscanf(line.c_str(), "b = %lf %lf", &Lx, &Ly);
			if(res != 3) {
				error_message << "Malformed headers found in an input configuration.\"b = <float> <float> <float>\" was expected, but \"" << line << "\" was found instead.";
				malformed_headers = true;
			}
		}
		if(malformed_headers) {
			double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
			res = sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, &t10, &t11, &t12, &t13, &t14, &t15);
			if(res == 15) {
				error_message << "\nSince the line contains 15 floats, likely the configuration contains more particles than specified in the topology file, or the header has been trimmed away.";
			}
			throw RCexception(error_message.str().c_str());
		}

		//std::getline(conf_input, line);
	}

	box->init(Lx, Ly);

	// the following part is always carried out in double precision since we want to be able to restart from confs that have
	// large numbers in the conf file and use float precision later
	int k, i;
	i = 0;
	std::string line;
	int conf_initial_line=9;
	while(!conf_input.eof() && i < N()) {
		BaseField *p = fields[i];
		std::getline(conf_input, line);
		auto spl_line = Utils::split_to_numbers(line, " ");

		p->init((int)spl_line[0], (int)spl_line[1]);
		p->CoM = std::vector<number> {spl_line[2], spl_line[3]};
		p->set_positions((int)spl_line[4], (int)spl_line[5], (int)spl_line[6]);

		p->nemQ = {spl_line[7], spl_line[8]};
		p->nemQ_old = {p->nemQ[0], p->nemQ[1]};
		p->Q00 = 0.5 * (p->nemQ[0] * p->nemQ[0] - p->nemQ[1] * p->nemQ[1]);
		p->Q01 = p->nemQ[0] * p->nemQ[1];

		for(k=0; k<p->subSize; k+=1){
			p->fieldScalar[k]=spl_line[conf_initial_line+(2*k)+1];
			p->area+=spl_line[conf_initial_line+(2*k)+1]*spl_line[conf_initial_line+(2*k)+1];
			p->sumF+=spl_line[conf_initial_line+(2*k)+1];
		}

		i++;
	}


	if(i != N()) {
		if(confs_to_skip > 0) {
			throw RCexception("Wrong number of particles (%d) found in configuration. Maybe you skipped too many configurations?", i);
		}
		else {
			throw RCexception("The number of lines found in configuration file (%d) doesn't match the parsed number of particles (%d)", i, N());
		}
	}

	interaction->check_input_sanity(fields);

	return true;
}


void SimBackend::apply_simulation_data_changes() {

}

void SimBackend::apply_changes_to_simulation_data() {

}

void SimBackend::add_output(ObservableOutputPtr new_output) {
	obs_outputs.push_back(new_output);
}

void SimBackend::remove_output(std::string output_file) {
	auto search = std::find_if(obs_outputs.begin(), obs_outputs.end(), [output_file](ObservableOutputPtr p) { return p->get_output_name() == output_file; });
	if(search == obs_outputs.end()) {
		throw RCexception("The output '%s' does not exist, can't remove it", output_file.c_str());
	}

	obs_outputs.erase(search);
}

void SimBackend::print_observables() {
	bool someone_ready = false;
	for(auto const &element : obs_outputs) {
		if(element->is_ready(current_step()))
			someone_ready = true;
	}

	if(someone_ready) {
		apply_simulation_data_changes();

		llint total_bytes = 0;
		for(auto const &element : obs_outputs) {
			if(element->is_ready(current_step())) {
				element->print_output(current_step());
			}
			total_bytes += element->get_bytes_written();
		}

		// here we control the timings; we leave the code a 30-second grace time to print the initial configuration
		double time_passed = (double) mytimer->get_time() / (double) CLOCKS_PER_SEC;
		if(time_passed > 30) {
			double MBps = (total_bytes / (1024. * 1024.)) / time_passed;
			OX_DEBUG("Current data production rate: %g MB/s (total bytes: %lld, time_passed: %g)" , MBps, total_bytes, time_passed);
			if(MBps > max_io) {
				throw RCexception(
						"Aborting because the program is generating too much data (%g MB/s).\n\t\t\tThe current limit is set to %g MB/s;\n\t\t\tyou can change it by setting max_io=<float> in the input file.",
						MBps, max_io);
			}
		}

		apply_changes_to_simulation_data();
	}

	backend_info = std::string("");
}

void SimBackend::update_observables_data() {
	bool updated = false;
	for(auto const &obs : config_info->observables) {
		if(obs->need_updating(current_step())) {
			if(!updated) {
				apply_simulation_data_changes();
				updated = true;
			}

			obs->update_data(current_step());
		}
	}
}


void SimBackend::print_conf(bool only_last) {
	apply_simulation_data_changes();

	if(!only_last) {
		obs_output_trajectory->print_output(current_step());
	}
	obs_output_last_conf->print_output(current_step());
}
