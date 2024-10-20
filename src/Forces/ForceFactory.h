#ifndef FORCEFACTORY_H_
#define FORCEFACTORY_H_

#include <string>

#include "../Utilities/RCexception.h"
#include "../Fields/BaseField.h"

#include "BaseForce.h"
#include "../Boxes/BaseBox.h"

/**
 * @brief Produces and adds external forces to the particles.
 *
 * This class is implemented as a singleton. See the comments in Logger.h/cpp for the singleton structure.
 *
 */

class ForceFactory {
private:
	static std::shared_ptr<ForceFactory> ForceFactoryPtr;
	ForceFactory();

public:
	virtual ~ForceFactory();

	static std::shared_ptr<ForceFactory> instance();

	/**
	 * @brief Produces and adds the force specified in the input file inp to the right particles.
	 *
	 * @param inp
	 * @param particles
	 * @param N
	 * @param box_side_ptr pointer to the box side. We use a pointer since the box size can change 
	 */
	void add_force(input_file &inp, std::vector<BaseField *> &fields, BaseBox *box_ptr);

	/// adds forces. Used by SimBackend and GeneratorManager
	void make_forces(std::vector<BaseField *> &fields, BaseBox *box_ptr);
};

#endif /* FORCEFACTORY_H_ */
