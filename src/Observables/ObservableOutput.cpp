#include <sstream>
#include <iostream>

#include "ObservableOutput.h"
#include "ObservableFactory.h"
#include "../Utilities/Utils.h"

using namespace std;

ObservableOutput::ObservableOutput(std::string &stream_string) :
				prefix("") {
	start_from = 0;
	stop_at = -1;
	bytes_written = 0;
	output_name = std::string("");
	input_file *obs_input = Utils::get_input_file_from_string(stream_string);
	input_file *sim_inp = CONFIG_INFO->sim_input;

	string out_name;
	getInputString(obs_input, "name", out_name, 1);
	// the prefix should be used only when the output is nor stdout nor stderr
	if(out_name != "stdout" && out_name != "stderr") {
		getInputString(sim_inp, "output_prefix", prefix, 0);
	}
	//sprintf(_output_name, "%s%s", _prefix.c_str(), out_name.c_str());
	output_name = prefix + out_name;

	getInputLLInt(obs_input, "start_from", &start_from, 0);
	getInputLLInt(obs_input, "stop_at", &stop_at, 0);

	int restart_step_counter = 0;
	getInputBoolAsInt(sim_inp, "restart_step_counter", &restart_step_counter, 0);
	append = !restart_step_counter;

	only_last = 0;
	getInputBoolAsInt(obs_input, "only_last", &only_last, 0);

	update_name_with_time = false;
	getInputBool(obs_input, "update_name_with_time", &update_name_with_time, 0);
	if(update_name_with_time) {
		base_name = output_name;
	}

	if(only_last && update_name_with_time) {
		throw oxDNAException("only_last and update_name_with_time are mutually exclusive");
	}

	linear = true;
	getInputBool(obs_input, "linear", &linear, 0);
	if(!linear) {
		getInputInt(obs_input, "log_ppc", &log_ppc, 1);
		getInputInt(obs_input, "log_n0", &log_n0, 1);
		getInputNumber(obs_input, "log_fact", &log_fact, 1);
		log_next = log_n0;
		log_pos_in_cycle = 1;
		log_n_cycle = 0;
		log_tot_cycle = (llint) round(log_n0 * pow(log_fact, (number) (log_ppc - 1.)));
	}
	else {
		getInputLLInt(obs_input, "print_every", &print_every, 1);
	}

	int i = 1;
	bool found = true;
	while(found) {
		stringstream ss;
		ss << "col_" << i;
		string obs_string;
		if(getInputString(obs_input, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			add_observable(obs_string);
		}
		else found = false;
		i++;
	}

	delete obs_input;
}

ObservableOutput::~ObservableOutput() {
	if(output_stream.is_open()) {
		output_stream.close();
	}
	clear();
}

void ObservableOutput::open_output() {
	if(output_stream.is_open()) {
		output_stream.close();
	}

	if(!strncmp(output_name.c_str(), "stderr", 512)) output = &std::cerr;
	else if(!strncmp(output_name.c_str(), "stdout", 512)) output = &std::cout;
	else {
		if(append && !update_name_with_time && !only_last) {
			output_stream.open(output_name.c_str(), ios_base::app);
		}
		else {
			output_stream.open(output_name.c_str());
		}

		output = &output_stream;
	}

	if(output->bad() || !output->good()) {
		throw oxDNAException("Stream %s not writable", output_name.c_str());
	}
}

void ObservableOutput::init() {
	if(!update_name_with_time) {
		open_output();
	}

	for(auto it = obss.begin(); it != obss.end(); it++) {
		(*it)->init();
	}
}

void ObservableOutput::clear() {
	obss.clear();
}

void ObservableOutput::add_observable(ObservablePtr new_obs) {
	obss.push_back(new_obs);
	CONFIG_INFO->observables.push_back(new_obs);
}

void ObservableOutput::add_observable(std::string obs_string) {
	std::shared_ptr<input_file> obs_inp(Utils::get_input_file_from_string(obs_string));

	ObservablePtr new_obs = ObservableFactory::make_observable(*obs_inp);
	add_observable(new_obs);
}

void ObservableOutput::change_output_file(string new_filename) {
	output_name = prefix + new_filename;
	open_output();
}

void ObservableOutput::set_next_log_step() {
	log_next = (llint) round((log_n0 * pow(log_fact, log_pos_in_cycle))) + log_tot_cycle * log_n_cycle;
	log_pos_in_cycle++;
	if(log_pos_in_cycle == log_ppc) {
		log_n_cycle++;
		log_pos_in_cycle = 0;
	}
}

bool ObservableOutput::is_ready(llint step) {
	if(stop_at > -1 && step > stop_at) {
		return false;
	}

	if(linear) {
		if(print_every < 1) return false;
		return ((step >= start_from) && (step % print_every == 0));
	}
	else {
		while(log_next < step) {
			set_next_log_step();
		}
		return (step == log_next);
	}
}

void ObservableOutput::print_output(llint step) {
	stringstream ss;
	for(auto it = obss.begin(); it != obss.end(); it++) {
		if(it != obss.begin()) ss << " ";
		ss << (*it)->get_output_string(step);
	}

	if(update_name_with_time) {
		string new_name = Utils::sformat("%s%lld", base_name.c_str(), step);
		change_output_file(new_name);
	}
	else if(only_last) open_output();

	ss << endl;
	std::string towrite = ss.str();
	bytes_written += (llint) towrite.length();
	*output << towrite;
	output->flush();

	if(only_last) {
		output_stream.close();
	}
	if(!linear) {
		set_next_log_step();
	}
}
