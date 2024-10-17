#ifndef BASEFIELD_H_
#define BASEFIELD_H_

#include <cstring>
#include <cstdlib>
#include <cassert>

#include <stdio.h>

#include "../defs.h"
#include "../Boxes/BaseBox.h"
#include "../Forces/BaseForce.h"

/**
 * @brief Base particle class. All particles must inherit from this class.
 */
class BaseField {
protected:
	void check();

public:
	std::vector<BaseForce *> ext_forces;

	BaseField();
	virtual ~BaseField();

	virtual void set_positions_initial(BaseBox *box) {}
	virtual void set_positions(BaseBox *box) {}
	virtual void set_positions(int offsetx, int offsety, int corner) {}
	virtual void copy_from(const BaseField &p);
	virtual void get_interaction_values(int R);

	virtual void init();
	virtual void init(int Lx, int Ly) {}

	std::vector<number> pos_shift = std::vector<number> {0., 0.};

	virtual int get_index() const {
		return index;
	}

	/**
	 * @brief Add an external force.
	 *
	 * @param f
	 * @return true if the force was added, false otherwise
	 */
	bool add_ext_force(BaseForce *f);
	std::vector<number> sf;

	inline void set_initial_forces(llint step, const std::shared_ptr<BaseBox> &box) {
		force = std::vector<number> {0., 0.};

		if(ext_forces.size() > 0) {
			std::vector<number> abs_pos = box->get_abs_pos(this);
			sf = std::vector<number> {0., 0.};
			for(auto ext_force : ext_forces) {
				sf = ext_force->value(step, abs_pos);
				force[0] += sf[0]; force[1] += sf[1];
			}
		}
	}

        /**
         * @brief Computes the interaction resulting from all the external forces acting on the particle. Stores the result in the ext_potential member.
         *
         * @param step current time step. Useful for forces that depend on time.
         * @param box pointer to the box object
         * @return true if the external potential was added, false otherwise
         */
        /// Total potential energy due to external forces
        double ext_potential;
        inline void set_ext_potential(llint step, BaseBox *box) {
                if(ext_forces.size() > 0) {
                        std::vector<number> abs_pos = box->get_abs_pos(this);
                        ext_potential = (double) 0.;
                        for(auto ext_force : ext_forces) {
                                ext_potential += ext_force->potential(step, abs_pos);
                        }
                }
        }


	/// Index of the particle. Usually it is a useful way of accessing arrays of particles
	int index;

	int type;

	/// External force exerted on the field
	std::vector<number> force;

	/// Positions of all interaction centers. This array must be initialized by child classes
	std::vector<number> CoM;

	/// Velocity of the particle
	number velocityX;
	number velocityY;
	std::vector<number> fieldDX;
	std::vector<number> fieldDY;

	number S01;
        number S00;
	std::vector<number> Fpressure;
	std::vector<number> Fshape;

	//child functions need to initialize everything below
	virtual int GetSubIndex(int site, BaseBox *box) = 0;
        virtual int GetSubXIndex(int site, BaseBox *box) = 0;
        virtual int GetSubYIndex(int site, BaseBox *box) = 0;
        std::vector<number> fieldScalar;
        std::vector<number> freeEnergy;
        std::vector<number> fieldScalar_old;
        std::vector<number> dfield_old;
	std::vector<int> neighbors_sub;

	int LsubX;
        int LsubY;
	int subSize;
	std::vector<int> offset;
        int sub_corner_bottom_left;

        number area;
        number sumF;
};

using FieldPtr = std::shared_ptr<BaseField>;

#endif /* BASEFIELD_H_ */
