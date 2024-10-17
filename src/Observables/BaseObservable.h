#ifndef BASEOBSERVABLE_H_
#define BASEOBSERVABLE_H_

#include "../defs.h"
#include "../Interactions/BaseInteraction.h"
#include "../Utilities/ConfigInfo.h"

/**
 * @brief Abstract class defining the basic interface for observables.
 */

class BaseObservable {
protected:
	/// Stores all the backend's information that may be needed by the observable
	ConfigInfo *config_info;

	std::string id;
	long long int_update_every = 0;
public:
	BaseObservable();

	virtual ~BaseObservable();

	virtual bool need_updating(llint curr_step);

	virtual void update_data(llint curr_step);

	/**
	 * @brief Returns the string to be printed in the output stream.
	 *
	 * @param curr_step
	 * @return observable output
	 */
	virtual std::string get_output_string(llint curr_step) = 0;

	/**
	 * @brief Reads the options stored in the input file associated to this observable and in the general input file.
	 *
	 * @param my_inp
	 * @param sim_inp
	 */
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	/**
	 * @brief Initializes the observable.
	 */
	virtual void init();

	void set_id(std::string idx) {
		id = idx;
	}

	std::string get_id() {
		return id;
	}
};

using ObservablePtr = std::shared_ptr<BaseObservable>;

#endif /* BASEOBSERVABLE_H_ */
