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
	virtual void update_positions(BaseBox *box) {}
	virtual void set_positions(int offsetx, int offsety, int corner, int corner_x, int corner_y, int size_x) {}
	virtual void set_properties_to_zero() {};
	virtual void copy_from(const BaseField &p);
	virtual void get_interaction_values(int R);
	virtual void check_borders(int q) {};
	//virtual void updateNeighborsSubSquarePeriodic(){return};

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
	number F_ext = 0;
	inline number set_F_ext(int q, number walls, number laplacian_walls) {
		if(ext_forces.size() > 0) {
			F_ext = 0.;
			for(auto ext_force : ext_forces) {
				F_ext += ext_force->free_energy(fieldScalar[q], walls, laplacian_walls);
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
        number ext_potential = 0;
        inline void set_ext_potential(int k, number walls) {
                if(ext_forces.size() > 0) {
                        ext_potential = 0.;
                        for(auto ext_force : ext_forces) {
                                ext_potential += ext_force->potential(k, walls);
                        }
                }
        }


	//makes changes after equilibration to the external forces of each phase field
        inline void set_ext_properties() {
                if(ext_forces.size() > 0) {
                        for(auto ext_force : ext_forces) {
                                ext_force->apply_changes_after_equilibration();
                        }
                }
        }


	//makes changes after equilibration to the external forces of each phase field
	std::vector<number> V_ext = std::vector<number> {0.,0.};
        inline void set_V_ext(int k, int lx, int ly) {
                if(ext_forces.size() > 0) {
                        for(auto ext_force : ext_forces) {
                                V_ext = ext_force->velocity_profile(k, lx, ly);
                        }
                }
        }


	/// Index of the particle. Usually it is a useful way of accessing arrays of particles
	int index;
	int type;
	int clock;

	/// Positions of all interaction centers. This array must be initialized by child classes
	std::vector<number> CoM;
	std::vector<number> CoM_old;
	//std::vector<number> tracer_particle;
	//std::vector<number> tracer_particle_old;

	/// Velocity of the particle
	std::vector<number> velocityX;
	std::vector<number> velocityY;
	std::vector<number> velocityX_correction;
	std::vector<number> phi_correction;
	number velocityX_CoM;
	number phi_CoM;
	//derivatives of the field
	std::vector<number> fieldDX;
	std::vector<number> fieldDY;
	std::vector<number> laplacianPhi;
	//derivatives of the field for the shape free energy
	std::vector<number> Phi00;
	std::vector<number> Phi01;

	//anisotropic terms
	std::vector<number> aniTerm1x;
	std::vector<number> aniTerm2x;
	std::vector<number> aniTerm1y;
	std::vector<number> aniTerm2y;


	//shape tensor
	number S01;
        number S00;

	//nematic tensor
	number Q00;
	number Q01;
	number NQ00;
	number NQ01;
	number thetaQ;
	number thetaQ_old;
	std::vector<number> nemQ;
	std::vector<number> nemQ_old;

	//passive and active forces
	std::vector<number> Fpassive_x;
	std::vector<number> Factive_x;
	std::vector<number> Fpassive_y;
	std::vector<number> Factive_y;
	std::vector<number> Pressure;

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

	//field patch values
	int LsubX;
        int LsubY;
	int subSize;
	std::vector<int> offset;
	std::vector<int> shear_velocity_sign;
        int sub_corner_bottom_left;
        number unrap_sub_corner_bottom_left_x;
        number unrap_sub_corner_bottom_left_y;
        int sub_corner_bottom_left_old;
	int border; //usually equals 4
	virtual void set_sub_border(){border=-1;};
	int x_sub_left, y_sub_top, x_sub_right, y_sub_bottom;
	std::vector<int> map_sub_to_box;
	std::vector<int> map_sub_to_box_x;
	std::vector<int> map_sub_to_box_y;

	std::vector<number> cos_x_table;
	std::vector<number> cos_y_table;
	std::vector<number> sin_x_table;
	std::vector<number> sin_y_table;

	//general properties of fields
        number area;
        number area_old;
        number perimeter;
        number sumF;
        number sumF_old;
};

using FieldPtr = std::shared_ptr<BaseField>;

#endif /* BASEFIELD_H_ */
