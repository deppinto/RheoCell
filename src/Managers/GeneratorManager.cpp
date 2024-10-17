#include "GeneratorManager.h"
#include "../Interactions/InteractionFactory.h"
#include "../Forces/ForceFactory.h"
#include "../Boxes/BoxFactory.h"

GeneratorManager::GeneratorManager(input_file inp, char *third_argument) :
				input(inp) {
	use_density = false;
	density = -1.;
	box_side = -1.;
	box_side_x = box_side_y;

	// first we look for an x in argv[2]
	std::string tmpstr(third_argument);
	size_t found_x = tmpstr.find('x');
	if(found_x != std::string::npos) {
		std::vector<std::string> sides_str = Utils::split(tmpstr, 'x');
		box_side_x = atof(sides_str[0].c_str());
		box_side_y = atof(sides_str[1].c_str());

		box_side = box_side_x;
		if(box_side_y < box_side) box_side = box_side_y;

		OX_LOG(Logger::LOG_INFO, "Detected non square box; box sides %g, %g", box_side_x, box_side_y);
		use_density = false;
	}
	else {
		number argv2 = atof(tmpstr.c_str());
		if(argv2 <= 0.) {
			throw oxDNAException("Refusing to run with box_side/density = %g (converted from \"%s\")", argv2, tmpstr.c_str());
		}
		if(argv2 < 2.) {
			use_density = true;
			density = argv2;
			OX_LOG(Logger::LOG_INFO, "Generating configuration with density %g", density);
		}
		else {
			use_density = false;
			box_side = argv2;
			box_side_x = box_side_y = box_side;
			OX_LOG(Logger::LOG_INFO, "Generating configuration with side %g", box_side);
		}
	}

	N = 0;
	interaction = NULL;

	external_forces = false;
	external_filename = std::string("");

	ConfigInfo::init(&fields);
}

GeneratorManager::~GeneratorManager() {
	for(auto p : fields) {
		delete p;
	}
}

void GeneratorManager::load_options() {
	CONFIG_INFO->sim_input = &input;

	getInputString(&input, "trajectory_file", trajectory, 1);
	getInputString(&input, "conf_file", output_conf, 1);

	// seed;
	int seed;
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
	OX_LOG(Logger::LOG_INFO, "Setting the random number generator with seed = %d", seed);
	srand48((long int) seed);

	Logger::instance()->get_settings(input);

	mybox = std::shared_ptr<BaseBox>(BoxFactory::make_box(input));

	interaction = InteractionFactory::make_interaction(input);
	interaction->get_settings(input);

	getInputBool(&input, "external_forces", &external_forces, 0);
}

void GeneratorManager::init() {
	interaction->init();

	N = interaction->get_N_from_topology();
	fields.resize(N);
	interaction->read_topology(fields);

	if(use_density) {
		box_side = pow(N / density, 1. / 2.);
		box_side_x = box_side_y = box_side;
		OX_LOG(Logger::LOG_INFO, "Generating square configuration with box_side %g", box_side);
	}
	else {
		density = N / (box_side_x * box_side_y);
		OX_LOG(Logger::LOG_INFO, "Generating configuration with density %g (%d particles, box sides %g %g)", density, N, box_side_x, box_side_y);
	}

	// setting the box of the interaction
	mybox->init(box_side_x, box_side_y);
	interaction->set_box(mybox.get());

	CONFIG_INFO->interaction = interaction.get();
	CONFIG_INFO->box = mybox.get();

	// initialise external forces
	ForceFactory::instance()->make_forces(fields, mybox.get());
}

void GeneratorManager::generate() {
	interaction->generate_random_configuration(fields);

	for(auto p : fields) {
		p->init();
		p->set_positions(CONFIG_INFO->box);
	} 

	std::ofstream conf_output(output_conf);
	conf_output.precision(15);

	conf_output << "t = 0" << std::endl;
	conf_output << "b = " << box_side_x << " " << box_side_y << std::endl;

	for(int i = 0; i < N; i++) {
		BaseField *p = fields[i];
		conf_output << p->LsubX << " "<< p->LsubY<< " "<< p->CoM[0] << " " << p->CoM[1] << " " << p->offset[0] << " " << p->offset[1] << " " << p->sub_corner_bottom_left << " ";
	        for(int q=0; q<p->subSize; q++)conf_output << p->GetSubIndex(q, CONFIG_INFO->box) << " " << p->fieldScalar[q] << " ";
		conf_output << std::endl;
	}
	conf_output.close();
}
