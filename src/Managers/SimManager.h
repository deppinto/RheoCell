#ifndef SIMMANAGER_H_
#define SIMMANAGER_H_

#include <cstring>
#include <ctime>

#include "../defs.h"
#include "../Backends/SimBackend.h"
#include "../Utilities/time_scales/time_scales.h"

/**
 * @brief Manages a simulation, be it on GPU or on CPU.
 *
 */
class SimManager {
protected:
	std::shared_ptr<SimBackend> backend;
	input_file input;
	time_scale time_scale_manager;
	int time_scale_var;
	llint steps, equilibration_steps;
	llint steps_run;
	int seed;
	int print_energy_every;
	int pid;
	int print_input;

public:
	SimManager(input_file inp);
	virtual ~SimManager();

	static bool stop;
	static bool started;
	virtual void load_options();
	virtual void init();
	virtual void run();
	virtual void analyse();
};

void gbl_termhandler(int);

#endif /* SIMMANAGER_H_ */
