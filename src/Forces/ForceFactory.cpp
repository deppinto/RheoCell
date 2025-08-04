#include "ForceFactory.h"

#include "ChannelWalls.h"
#include "ShearFlowChannel.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <sstream>

using namespace std;

std::shared_ptr<ForceFactory> ForceFactory::ForceFactoryPtr = nullptr;

ForceFactory::ForceFactory() {

}

ForceFactory::~ForceFactory() {

}

std::shared_ptr<ForceFactory> ForceFactory::instance() {
	// we can't use std::make_shared because ForceFactory's constructor is private
	if(ForceFactoryPtr == nullptr) {
		ForceFactoryPtr = std::shared_ptr<ForceFactory>(new ForceFactory());
	}

	return ForceFactoryPtr;
}

void ForceFactory::add_force(input_file &inp, std::vector<BaseField *> &fields, BaseBox * box_ptr) {
	string type_str;
	getInputString(&inp, "type", type_str, 1);

	ForcePtr extF;

	if(type_str.compare("channel_walls") == 0) extF = std::make_shared<ChannelWalls>();
	else if(type_str.compare("moving_walls") == 0) extF = std::make_shared<ChannelWalls>();
	else if(type_str.compare("circle_walls") == 0) extF = std::make_shared<ChannelWalls>();
	else if(type_str.compare("triangle_walls") == 0) extF = std::make_shared<ChannelWalls>();
	else if(type_str.compare("shear_flow_channel") == 0) extF = std::make_shared<ShearFlowChannel>();
	else throw RCexception("Invalid force type `%s\'", type_str.c_str());

	string group = string("default");
	getInputString(&inp, "group_name", group, 0);
	extF->set_group_name(group);

	std::vector<int> field_ids;
	std::string description;
	std::tie(field_ids, description) = extF->init(inp);

	CONFIG_INFO->add_force_to_fields(extF, field_ids, description);
}

void ForceFactory::make_forces(std::vector<BaseField *> &fields, BaseBox *box) {
	bool external_forces = false;
	getInputBool(CONFIG_INFO->sim_input, "external_forces", &external_forces, 0);

	if(external_forces) {
		std::string external_filename;
		getInputString(CONFIG_INFO->sim_input, "external_forces_file", external_filename, 1);

		ifstream external(external_filename.c_str());
		if(!external.good ()) {
			throw RCexception ("Can't read external_forces_file '%s'", external_filename.c_str());
		}

		bool is_json = false;
		getInputBool(CONFIG_INFO->sim_input, "external_forces_as_JSON", &is_json, 0);

		if(is_json) {
			OX_LOG(Logger::LOG_INFO, "Parsing JSON force file %s", external_filename.c_str());

			nlohmann::json my_json = nlohmann::json::parse(external);

			for(auto &force_json : my_json) {
				input_file force_input;
				force_input.init_from_json(force_json);
				ForceFactory::instance()->add_force(force_input, fields, box);
			}
		}
		else {
			OX_LOG(Logger::LOG_INFO, "Parsing force file %s", external_filename.c_str());

			//char line[512], typestr[512];
			int open, justopen, a;

			justopen = open = 0;
			a = external.get();
			bool is_commented = false;
			stringstream external_string;
			while(external.good()) {
				justopen = 0;
				switch(a) {
					case '#':
					is_commented = true;
					break;
					case '\n':
					is_commented = false;
					break;
					case '{':
					if(!is_commented) {
						open++;
						justopen = 1;
					}
					break;
					case '}':
					if(!is_commented) {
						if(justopen) throw RCexception ("Syntax error in '%s': nothing between parentheses", external_filename.c_str());
						open--;
					}
					break;
					default:
					break;
				}

				if(!is_commented) external_string << (char)a;
				if(open > 1 || open < 0) throw RCexception ("Syntax error in '%s': parentheses do not match", external_filename.c_str());
				a = external.get();
			}
			external_string.clear();
			external_string.seekg(0, ios::beg);
			a = external_string.get();
			while(external_string.good()) {
				while (a != '{' && external_string.good()) {
					a = external_string.get();
				}
				if(!external_string.good()) {
					break;
				}

				a = external_string.get();
				std::string input_string("");
				while (a != '}' && external_string.good()) {
					input_string += a;
					a = external_string.get();
				}
				input_file input;
				input.init_from_string(input_string);

				ForceFactory::instance()->add_force(input, fields, box);
			}
		}

		external.close();
		OX_LOG(Logger::LOG_INFO, "   Force file parsed");
	}
}
