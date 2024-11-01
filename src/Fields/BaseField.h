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
	virtual void set_properties_to_zero() {};
	virtual void copy_from(const BaseField &p);
	virtual void get_interaction_values(int R);
	virtual void check_borders(int q, int box_size_x, int box_size_y) {};

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
	number F_ext;
	inline number set_F_ext(int k, number phi, number walls) {
		if(ext_forces.size() > 0) {
			F_ext = 0.;
			for(auto ext_force : ext_forces) {
				F_ext += ext_force->free_energy(k, phi, walls);
			}
		}
		return F_ext;
	}

        /**
         * @brief Computes the interaction resulting from all the external forces acting on the particle. Stores the result in the ext_potential member.
         *
         * @param step current time step. Useful for forces that depend on time.
         * @param box pointer to the box object
         * @return true if the external potential was added, false otherwise
         */
        /// Total potential energy due to external forces
        number ext_potential;
        inline void set_ext_potential(int k, number walls) {
                if(ext_forces.size() > 0) {
                        ext_potential = 0.;
                        for(auto ext_force : ext_forces) {
                                ext_potential += ext_force->potential(k, walls);
                        }
                }
        }


	/// Index of the particle. Usually it is a useful way of accessing arrays of particles
	int index;
	int type;

	/// Positions of all interaction centers. This array must be initialized by child classes
	std::vector<number> CoM;

	/// Velocity of the particle
	std::vector<number> velocityX;
	std::vector<number> velocityY;
	//derivatives of the field
	std::vector<number> fieldDX;
	std::vector<number> fieldDY;

	//shape tensor
	number S01;
        number S00;

	//nematic tensor
	number Q00;
	number Q01;
	number thetaQ;
	number thetaQ_old;
	std::vector<number> nemQ;
	std::vector<number> nemQ_old;

	//passive and active forces
	std::vector<number> Fpassive;
	std::vector<number> Factive;

	//child functions need to initialize everything below
	virtual int GetSubIndex(int site, BaseBox *box) = 0;
        virtual int GetSubXIndex(int site, BaseBox *box) = 0;
        virtual int GetSubYIndex(int site, BaseBox *box) = 0;

	//vectors for integration
        std::vector<number> fieldScalar;
        std::vector<number> freeEnergy;
        std::vector<number> fieldScalar_old;
        std::vector<number> dfield_old;
	std::vector<int> neighbors_sub;

	//fiedl patch values
	int LsubX;
        int LsubY;
	int subSize;
	std::vector<int> offset;
        int sub_corner_bottom_left;
	int border;
	int x_sub_left, y_sub_top, x_sub_right, y_sub_bottom;
	std::vector<int> map_sub_to_box;
	std::vector<int> map_sub_to_box_x;
	std::vector<int> map_sub_to_box_y;

	//general properties of fields
        number area;
        number sumF;
};

using FieldPtr = std::shared_ptr<BaseField>;

#endif /* BASEFIELD_H_ */
