#ifndef MULTIPHASEFIELD_H_
#define MULTIPHASEFIELD_H_

#include "BaseField.h"
#include "../Boxes/BaseBox.h"

/**
 * @brief A customisable particle. Used by CustomInteraction.
 */

class MultiPhaseField: public BaseField {
protected:
	int init_radius=4;
	int init_radius2=4*4;
	void setNeighborsSubSquareDirichlet();
	void setNeighborsSubSquarePeriodic();
	void resizing();

public:
	MultiPhaseField();
	MultiPhaseField(int Lx, int Ly);
	virtual ~MultiPhaseField();

	void setNeighborsSub();

	void set_positions_initial(BaseBox *box);
	void set_positions(BaseBox *box);
	void set_positions(int offsetx, int offsety, int corner);
	void set_properties_to_zero();
	void init();
	virtual void init(int Lx, int Ly);
	void get_interaction_values(int R); 
	void check_borders(int q, int box_size_x, int box_size_y);

	int GetSubIndex(int site, BaseBox *box);
	int GetSubXIndex(int site, BaseBox *box);
	int GetSubYIndex(int site, BaseBox *box);
};

#endif /* MULTIPHASEFIELD_H_ */
