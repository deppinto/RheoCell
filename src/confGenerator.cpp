#include "defs.h"
#include "Managers/GeneratorManager.h"
#include "Utilities/SignalManager.h"

void print_version() {
	exit(-1);
}

int main(int argc, char *argv[]) {
	try {
		Logger::init();
		SignalManager::manage_segfault();
		if(argc < 3) {
			throw oxDNAException("Usage is '%s input_file [box_size|density]'\nthe third argument will be interpreted as a density if it is less than 2.0", argv[0]);
		}
		else if(argc > 1 && !strcmp(argv[1], "-v")) {
			print_version();
		}

		input_file input(true);
		input.init_from_command_line_args(argc, argv, 1);

		GeneratorManager mygenerator(input, argv[2]);
		OX_DEBUG("Loading options");
		mygenerator.load_options();

		OX_DEBUG("Initializing");
		mygenerator.init();

		OX_DEBUG("Running");
		OX_LOG(Logger::LOG_INFO, "Starting generation");
		mygenerator.generate();

		OX_LOG(Logger::LOG_INFO, "END OF THE GENERATION, everything went OK!");
	}
	catch (oxDNAException &e) {
		OX_LOG(Logger::LOG_ERROR, "%s", e.what());
		exit(1);
	}

	return 0;
}
