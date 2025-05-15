/*
 * BaseInteraction is a class that is meant "contain" all the interaction
 * types
 */

#ifndef BASE_INTERACTION_H
#define BASE_INTERACTION_H

#include "../defs.h"
#include "../Fields/BaseField.h"
#include "../Boxes/BaseBox.h"
#include "../Utilities/Utils.h"

#include <fstream>
#include <set>
#include <vector>
#include <stdio.h>
#include <string.h>

/**
 * @brief Base class for managing particle-particle interactions. It is an abstract class.
 */
class BaseInteraction {
private:

protected:
	BaseBox *box;

	int R;
	number friction = 1.0;
        number rcut, sqr_rcut;
	/// This is useful for "hard" potentials
	bool is_infinite;
	
	number computed_r;
	char topology_filename[256];

	number energy_threshold;

        //!Allow a openMP threads
        int omp_thread_num = 1;
        //set number of threads
        virtual void set_omp_threads(int num_threads){omp_thread_num = num_threads;};
        virtual void set_omp_tasks(int num_threads){};

public:
	/**
	 * @brief Basic constructor. By default, it does not need anything.
	 */
	BaseInteraction();
	virtual ~BaseInteraction();

	virtual void set_box(BaseBox *boxArg) {
		box = boxArg;
	}

	//virtual void get_settings(input_file &inp);

	/**
	 * Initialization of class constants.
	 */
	virtual void init() = 0;

	virtual void apply_changes_after_equilibration(){};


	/**
	 * @brief Handles particle allocation. Child classes must implement it.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void allocate_fields(std::vector<BaseField *> &fields) = 0;

	/**
	 * @brief Check whether the initial configuration makes sense.
	 *
	 * Since this operation is interaction-dependent, each interaction must provide its own implementation.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void check_input_sanity(std::vector<BaseField *> &fields) = 0;

	/**
	 * @brief Computes the total interaction between particles p and q.
	 *
	 * If r is not given or NULL, it is computed from scratch. It can optionally update forces and torques exerted on p and q.
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return pair-interaction energy
	 */
	//virtual number field_interaction(BaseField *p, int q, bool update_forces = false) = 0;

	/**
	 * @brief Returns the total potential energy of the system
	 *
	 * @param particles
	 * @param N
	 * @return
	 */
	virtual number get_system_energy(std::vector<BaseField *> &fields);

	/**
	 * @brief Returns the state of the interaction
	 */
	bool get_is_infinite() {
		return is_infinite;
	}

	/**
	 * @brief Sets _is_infinite
	 *
	 * @param arg bool
	 */
	void set_is_infinite(bool arg) {
		is_infinite = arg;
	}

	/**
	 * @brief Check whether the two particles overlaps. Used only in the generation of configurations. Can be overloaded.
	 *
	 * The default implementation does not allow two particles to have an
	 * interaction strength greater than 100. More subtle implementations
	 * may be needed in special cases.
	 *
	 * @param p
	 * @param q
	 * @return bool whether the two particles overlap
	 */
	virtual bool generate_random_configuration_overlap(BaseField *p, BaseField *q);

	/**
	 * @brief Generate an initial configuration. Can be overloaded.
	 *
	 * The default function creates a random configuration in the most simple way possible:
	 * puts particles in one at a time, and using {\@link generate_random_configuration_overlap}
	 * to test if two particles overlap.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void generate_random_configuration(std::vector<BaseField *> &fields);
	virtual void generate_lattice_configuration(std::vector<BaseField *> &fields);
	virtual void generate_cluster_configuration(std::vector<BaseField *> &fields);


	/**
         * @brief Signals the interaction that an energy computation is about to begin.
         *
         * By default this method does nothing, but interactions inheriting from this interface may
         * need to initialise or reset other data structures before computing the energy or the force acting on all particles.
         *
         */
        virtual void begin_energy_computation();
	virtual void begin_energy_computation(std::vector<BaseField *> &fields);
	virtual void resetSums(int k);
	virtual void updateFieldProperties(BaseField *p, int q, int k);
	virtual void updateDirectedActiveForces(number dt, BaseField *p, bool store){};


        /**
         * @brief Read the topology of the interactions. The defaut
         * method is empty.
         *
         * This function is set up so that the topology of the interactions
         * are read in the code. The different interactions might require
         * different formats, so it is implemented here. The default
         *
         * @param N_strands Pointer to an int which will be filled with the
         * number of "strands" (unbreakable clusters of connected particles)
         * @param particles array of particles.
         */
        virtual void read_topology(std::vector<BaseField *> &fields);

        /**
         * @brief Returns the number of particles, as written in the topology.
         *
         * @return number of particles
         */
        virtual int get_N_from_topology();

	virtual void update_sub_to_box_map(BaseField *p, int q, int sub_site, int sub_site_x, int sub_site_y);
	virtual number get_velocity_x(BaseField *p, int q){return p->velocityX[q];}
	virtual number get_velocity_y(BaseField *p, int q){return p->velocityY[q];}
	virtual number get_total_force_x(BaseField *p, int q){return p->velocityX[q] * friction;}
	virtual number get_total_force_y(BaseField *p, int q){return p->velocityY[q] * friction;}

	virtual void get_settings(input_file &inp);
	number K;
	number U;
};

using InteractionPtr = std::shared_ptr<BaseInteraction>;

#endif /* BASE_INTERACTION_H */
