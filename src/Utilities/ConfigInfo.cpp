#include "ConfigInfo.h"

#include "oxDNAException.h"
#include "../Fields/BaseField.h"
#include "../Interactions/BaseInteraction.h"
#include "../Forces/BaseForce.h"
#include "../Observables/BaseObservable.h"

std::shared_ptr<ConfigInfo> ConfigInfo::config_info = nullptr;

ConfigInfo::ConfigInfo(std::vector<BaseField *> *ps) :
				fields_pointer(ps)
				 {

}

ConfigInfo::~ConfigInfo() {

}

void ConfigInfo::update_temperature(number new_T) {
	temperature = new_T;
	notify("T_updated");
}

void ConfigInfo::add_force_to_fields(std::shared_ptr<BaseForce> force, std::vector<int> field_ids, std::string force_description) {
	forces.push_back(force);

	if(field_ids[0] != -1) {
		for(auto id : field_ids) {
			fields()[id]->add_ext_force(force.get());
			OX_LOG(Logger::LOG_INFO, "Adding a %s on particle %d", force_description.c_str(), id);
		}
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding a %s on ALL particles", force_description.c_str());
		for(auto p: fields()) {
			p->add_ext_force(force.get());
		}
	}
}

ForcePtr ConfigInfo::get_force_by_id(std::string id) {
	auto it = std::find_if(forces.begin(), forces.end(), [&id](const ForcePtr& obj) {return obj->get_id() == id;});

	if(it != forces.end()) {
		return *it;
	}

	return nullptr;
}

ObservablePtr ConfigInfo::get_observable_by_id(std::string id) {
	auto it = std::find_if(observables.begin(), observables.end(), [&id](const ObservablePtr& obj) {return obj->get_id() == id;});

	if(it != observables.end()) {
		return *it;
	}

	return nullptr;
}

void ConfigInfo::subscribe(std::string event, std::function<void()> callback) {
	event_callbacks[event].emplace_back(callback);
}

void ConfigInfo::notify(std::string event) {
	for(auto callback : event_callbacks[event]) {
		callback();
	}
}

void ConfigInfo::set(BaseInteraction *i, std::string *info, BaseBox *abox) {
	interaction = i;
	backend_info = info;
	box = abox;
}

void ConfigInfo::init(std::vector<BaseField *> *ps) {
	if(config_info != nullptr) {
		throw oxDNAException("The ConfigInfo object has been already initialised");
	}

	config_info = std::shared_ptr<ConfigInfo>(new ConfigInfo(ps));
}

void ConfigInfo::clear() {
	config_info.reset();
	config_info = nullptr;
}

const FlattenedConfigInfo &ConfigInfo::flattened_configuration() {
	flattened_conf.update(curr_step, fields());
	return flattened_conf;
}
