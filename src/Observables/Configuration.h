#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "BaseObservable.h"

/**
 * @brief Prints ascii configurations. It can be extended to provide further
 * output types (e.g. for visualisation).
 *
 */

class Configuration: public BaseObservable {
protected:
	int only_type;
	std::set<int> visible_fields;
        std::set<int> hidden_fields;

	/**
	 * @brief Returns the configuration header(s)
	 *
	 * @param step
	 * @return
	 */
	virtual std::string headers(llint step);

	/**
	 * @brief Returns the portion of output for the given particle
	 *
	 * @param p
	 * @return
	 */
	virtual std::string field(BaseField *p);

	/**
	 * @brief Returns the configuration output for the whole system. It does not comprise the headers.
	 *
	 * @param step
	 * @return
	 */
	virtual std::string configuration(llint step);

public:
	Configuration();
	virtual ~Configuration();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	std::string get_output_string(llint curr_step);
};

#endif /* CONFIGURATION_H_ */

