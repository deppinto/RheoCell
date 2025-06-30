#ifndef BASEFORCE_H_
#define BASEFORCE_H_

#include <string>
#include "../defs.h"
#include "../Utilities/RCexception.h"
#include "../Utilities/Utils.h"
#include "../Utilities/ConfigInfo.h"

// forward declarations of BaseParticle and BaseBox; needed to compile
class BaseField;
class BaseBox;

/**
 * @brief Base class for external forces. All external forces inherit from here.
 *
 * Note: This class contains public members with names starting with underscores.
 * We have to change scope policy due to the fact that GPU classes
 * require access to these members.
 *
 */

class BaseForce {
private:
	/**
	 * @brief A name of the group this force belongs to.
	 *
	 * Different forces can be grouped under the same name. This can be used to
	 * act separately on different forces. For example, it can be used together
	 * with the observable ForceEnergy to print only the energy due to specific
	 * groups of forces.
	 */
	std::string group_name = "default";
	std::string id = "";
	std::string type = "type_unread";

public:
	/**
	 * @brief standard members for forces
	 *
	 * we need these members to be public because
	 * we need access in order to copy these numbers
	 * to the GPU memory
	 */
	BaseField *p_ptr;

	BaseForce();
	virtual ~BaseForce();

	/** 
	 * @brief init function
	 *
	 * @return the list of particles it will act on and the force description
	 *
	 * This function initialises the force object and returns the list of particles it will act on
	 */
	virtual std::tuple<std::vector<int>, std::string> init(input_file &inp);

	void set_group_name(std::string &name) {
		group_name = name;
	}

	std::string get_group_name() {
		return group_name;
	}

	void set_id(std::string idx) {
		id = idx;
	}

	std::string get_id() {
		return id;
	}

	std::string get_type() {
		return type;
	}

	/**
	 * @brief returns value of the force (a vector)
	 *
	 * @param step useful for forces that depend on time
	 * @param pos position of the particle
	 */
	virtual number free_energy(number phi, number walls, number laplacianSquare){return 0;};

	/**
	 * @brief returns value of the potential associated to the force (a number)
	 *
	 * @param step useful for forces that depend on time
	 * @param pos position of the particle
	 */
	virtual number potential(int k, number walls){return 0;};

	virtual void apply_changes_after_equilibration() {};
	virtual std::vector<number> velocity_profile(int k, int lx, int ly) {return std::vector<number>{0.,0.};};
};

using ForcePtr = std::shared_ptr<BaseForce>;

#endif /* BASEFORCE_H_ */
