#ifndef SIMBACKEND_H_
#define SIMBACKEND_H_

#define SIM_PC 0

#include "../defs.h"
#include "../Observables/ObservableOutput.h"

#include <cmath>
#include <fstream>
#include <cfloat>
#include <vector>
#include <map>

class BaseInteraction;
class BaseBox;
class BaseField;
class ConfigInfo;
class Timer;

/**
 * @brief Base backend class. Every backend inherits from this class.
 *
 */

class SimBackend {
protected:
	std::string backend_info;

	std::shared_ptr<Timer> mytimer;

	number max_io;

	bool reseed;

	std::shared_ptr<BaseBox> box;

	int conf_interval;
	std::string conf_filename;
	bool custom_conf_name;
	char custom_conf_str[256];
	std::ifstream conf_input;
	llint read_conf_step;
	bool restart_step_counter;

	/// Vector of ObservableOutput used to manage the simulation output
	std::vector<ObservableOutputPtr> obs_outputs;
	ObservableOutputPtr obs_output_stdout;
	ObservableOutputPtr obs_output_file;
	ObservableOutputPtr obs_output_trajectory;
	ObservableOutputPtr obs_output_last_conf;
	ObservableOutputPtr obs_output_equilibrated_conf;

	/// Shared pointer to the interaction manager
	InteractionPtr interaction;

	/// potential energy
        number U;
	/// kinetic energy
        number K;


	/// array of pointers to particle objects
	std::vector<BaseField *> fields;

	/// object that stores pointers to a few important variables that need to be shared with other objects
	ConfigInfo *config_info = nullptr;

	int N_updates;
	int confs_to_skip;

public:
	SimBackend();
	virtual ~SimBackend();

	/**
	 * @brief Loads all the required settings from the input file.
	 *
	 * @param inp
	 */
	virtual void get_settings(input_file &inp);
	/**
	 * @brief Initializes the backend.
	 *
	 * @param path
	 */
	virtual void init();

	/**
	 * @brief Reads the next configuration from the conf_file.
	 *
	 * It can read it either in binary or ascii format.
	 *
	 * @param binary whether conf_file is to be parsed in ascii or binary format
	 * @return true if the operation was successful, false otherwise
	 */
	virtual bool read_next_configuration(bool binary=false);

	int N() {
		return fields.size();
	}

	/**
	 * @brief Returns how many times verlet lists have been updated.
	 *
	 * @return number of Verlet lists updates
	 */
	int get_N_updates() {
		return N_updates;
	}

	void add_output(ObservableOutputPtr new_output);
	void remove_output(std::string output_file);

	/**
	 * @brief Prints the observables attached to the backend.
	 */
	virtual void print_observables();

	virtual void update_observables_data();

	virtual void print_conf(bool only_last=false);
	virtual void print_equilibration_info(bool only_last=false);

	/**
	 * @brief Performs a simulation step.
	 */
	virtual void sim_step() = 0;

	/**
         * @brief Performs an equilibration step.
         */
        virtual void eq_step() = 0;

	/**
	 * @brief Synchronize the simulation data with the data structures that are used to analyse/print the current simulation status.
	 */
	virtual void apply_simulation_data_changes();

	/**
	 * @brief Update the simulation data, so that changes done to the data structures are taken into account by the simulation engine.
	 */
	virtual void apply_changes_to_simulation_data();

	long long int current_step() {
		return config_info->curr_step;
	}

	void increment_current_step() {
		config_info->curr_step++;
	}

	llint start_step_from_file;
};

#endif /* SIMBACKEND_H_ */
