#ifndef LEBCMULTIPHASEFIELD_H_
#define LEBCMULTIPHASEFIELD_H_

#include "BaseField.h"
#include "../Boxes/BaseBox.h"

/**
 * @brief A customisable particle. Used by CustomInteraction.
 */

class LEBcMultiPhaseField: public BaseField {
protected:
	int init_radius=4;
	int init_radius2=4*4;
	void setNeighborsSubSquareDirichlet();
	void setNeighborsSubSquarePeriodic();
	void resizing();

public:
	LEBcMultiPhaseField();
	LEBcMultiPhaseField(int Lx, int Ly);
	virtual ~LEBcMultiPhaseField();

	void setNeighborsSub();

	void set_positions_initial(BaseBox *box);
	void set_positions(BaseBox *box);
	void update_positions(BaseBox *box);
	void set_positions(int offsetx, int offsety, int corner, int corner_x, int corner_y, int size_x);
	void set_properties_to_zero();
	void init();
	virtual void init(int Lx, int Ly);
	void get_interaction_values(int R); 
	void check_borders(int q);

	int GetSubIndex(int site, BaseBox *box);
	int GetSubXIndex(int site, BaseBox *box);
	int GetSubYIndex(int site, BaseBox *box);

	virtual void set_sub_border(){border=6;};
};

#endif /* LEBCMULTIPHASEFIELD_H_ */
