#ifndef GENERATORMANAGER_H_
#define GENERATORMANAGER_H_

#include "../defs.h"
#include "../Backends/SimBackend.h"
#include "../Fields/BaseField.h"
#include "../Boxes/BaseBox.h"

/**
 * @brief Manages the generation of an initial configuration.
 */
class GeneratorManager {
protected:
	input_file input;
	char output_conf[256];
	char trajectory[256];

	InteractionPtr interaction;

	bool use_density;
	double box_side;
	double box_side_x, box_side_y;
	double density;
	int N;
	std::vector<BaseField *> fields;

	bool external_forces;
	std::string external_filename;

	std::shared_ptr<BaseBox> mybox;

public:
	GeneratorManager(input_file inp, char *third_argument);
	virtual ~GeneratorManager();

	void load_options();
	void init();
	void generate();
};

#endif /* GENERATORMANAGER_H_ */
