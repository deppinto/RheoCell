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

public:
	MultiPhaseField();
	MultiPhaseField(int Lx, int Ly);
	virtual ~MultiPhaseField();

	void setNeighborsSub();

	void set_positions_initial(BaseBox *box);
	void set_positions(BaseBox *box);
	void set_positions(int offsetx, int offsety, int corner);
	void init();
	virtual void init(int Lx, int Ly);
	void get_interaction_values(int R); 

	int GetSubIndex(int site, BaseBox *box);
	int GetSubXIndex(int site, BaseBox *box);
	int GetSubYIndex(int site, BaseBox *box);
};

#endif /* MULTIPHASEFIELD_H_ */
