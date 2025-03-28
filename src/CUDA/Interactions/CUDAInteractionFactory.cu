#include <string>

#include "CUDAInteractionFactory.h"

#include "CUDAWetModel.h"
#include "../../Utilities/Utils.h"

CUDAInteractionFactory::CUDAInteractionFactory() {

}

CUDAInteractionFactory::~CUDAInteractionFactory() {

}

std::shared_ptr<CUDABaseInteraction> CUDAInteractionFactory::make_interaction(input_file &inp) {
	// The default interaction is the wet model
	string inter_type("wetmodel");
	getInputString(&inp, "interaction_type", inter_type, 0);

	if(!inter_type.compare("wetmodel")) return std::make_shared<CUDAWetModel>();
	else {
		std::string cuda_name(inter_type);
		cuda_name = "CUDA" + cuda_name;
		std::shared_ptr<CUDABaseInteraction> res = std::dynamic_pointer_cast<CUDABaseInteraction>(PluginManager::instance()->get_interaction(cuda_name));
		if(res == nullptr) {
			throw oxDNAException ("CUDA interaction '%s' not found. Aborting", cuda_name.c_str());
		}
		return res;
	}
}
