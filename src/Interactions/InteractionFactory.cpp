#include "InteractionFactory.h"

#include "SimpleMultiField.h"
#include "ActiveMultiField.h"

InteractionPtr InteractionFactory::make_interaction(input_file &inp) {
	std::string inter_type("simplefield");
	getInputString(&inp, "interaction_type", inter_type, 0);
	std::string backend;
	getInputString(&inp, "backend", backend, 0);

	if(inter_type.compare("simplefield") == 0) {
		if(backend.compare("CUDA") == 0) return std::make_shared<SimpleMultiField>();
		else return std::make_shared<SimpleMultiField>();
	}
	else if(inter_type.compare("activefield") == 0){
                if(backend.compare("CUDA") == 0) return std::make_shared<ActiveMultiField>();
                else return std::make_shared<ActiveMultiField>();
        }

	else {
		throw RCexception("Interaction '%s' not found. Aborting", inter_type.c_str());
	}
}
