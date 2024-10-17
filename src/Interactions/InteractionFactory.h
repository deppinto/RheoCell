#ifndef INTERACTIONFACTORY_H_
#define INTERACTIONFACTORY_H_

#include "BaseInteraction.h"
#include "../Utilities/parse_input/parse_input.h"

/**
 * @brief Static factory class. Its only public method builds an {@link IBaseInteraction interaction}.
 *
 */
class InteractionFactory {
public:
	InteractionFactory() = delete;
	virtual ~InteractionFactory() = delete;

	/**
	 * @brief Builds the interaction instance.
	 *
	 * @param inp
	 * @return a pointer to the newly built interaction
	 */
	static InteractionPtr make_interaction(input_file &inp);
};

#endif /* INTERACTIONFACTORY_H_ */
