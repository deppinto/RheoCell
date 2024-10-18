#include "defs.h"
#include "Managers/SimManager.h"
#include "Utilities/SignalManager.h"
#include "Utilities/RCexception.h"
#include "Utilities/Timings.h"

using namespace std;

void print_version() {
	exit(-1);
}

int main(int argc, char *argv[]) {
	try {
		Logger::init();
		SignalManager::manage_segfault();
		TimingManager::init();

		if(argc < 2) {
			throw RCexception("Usage is '%s input_file'", argv[0]);
		}
		if(!strcmp(argv[1], "-v")) {
			print_version();
		}

		input_file input(true);
		input.init_from_command_line_args(argc, argv);

		SimManager mysim(input);
		mysim.load_options();

		OX_DEBUG("Initializing");
		mysim.init();

		OX_DEBUG("Running");
		mysim.run();

		OX_LOG(Logger::LOG_INFO, "END OF THE SIMULATION, everything went OK!");
	}
	catch (RCexception &e) {
		OX_LOG(Logger::LOG_ERROR, "%s", e.what());
		return 1;
	}

	TimingManager::clear();

	OX_LOG(Logger::LOG_NOTHING, "");
	return 0;
}

