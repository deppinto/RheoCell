#include "ObservableFactory.h"

#include "Configuration.h"

#include <nlohmann/json.hpp>


ObservablePtr ObservableFactory::make_observable(input_file &obs_inp) {
	char obs_type[512];
	getInputString(&obs_inp, "type", obs_type, 1);

	ObservablePtr res = nullptr;

	if(!strncasecmp(obs_type, "configuration", 512)) res = std::make_shared<Configuration>();
	else {
		if(res == NULL) throw RCexception("Observable '%s' not found. Aborting", obs_type);
	}

	res->get_settings(obs_inp, *(CONFIG_INFO->sim_input));

	return res;
}

std::vector<ObservableOutputPtr> ObservableFactory::make_observables(std::string prefix) {
	std::vector<ObservableOutputPtr> result;

	std::string obs_input_base = Utils::sformat("%sdata_output_", prefix.c_str());
	std::string obs_key = Utils::sformat("%sobservables_file", prefix.c_str());

	int i = 1;
	bool found = true;
	while(found) {
		std::stringstream ss;
		ss << obs_input_base << i;
		std::string obs_string;
		if(getInputString(CONFIG_INFO->sim_input, ss.str().c_str(), obs_string, 0) == KEY_FOUND) {
			auto new_obs_out = std::make_shared<ObservableOutput>(obs_string);
			result.push_back(new_obs_out);
		}
		else {
			found = false;
		}

		i++;
	}

	return result;
}
