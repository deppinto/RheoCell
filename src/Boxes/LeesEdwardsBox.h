#ifndef LEESEDWARDSBOX_H_
#define LEESEDWARDSBOX_H_

#include "BaseBox.h"

class LeesEdwardsBox: public BaseBox {
protected:
	number factor;
        int sidex, sidey;
        std::vector<int> sides;
	llint last_step = -1;
	number delta_x = 0.;

public:
	LeesEdwardsBox();
	virtual ~LeesEdwardsBox();

	virtual void get_settings(input_file &inp);
	virtual void init(int Lx, int Ly);

	virtual int getElementX(int site, int distX);
        virtual int getElementY(int site, int distY);
        virtual int getElement(int site, int distX, int distY);

        virtual int getXsize();
        virtual int getYsize();
        virtual int getBoxsize();

	virtual std::vector<number> min_image(const std::vector<number> &v1, const std::vector<number> &v2) const;
	virtual std::vector<number> min_image_PBC(const std::vector<number> &v1, const std::vector<number> &v2) const;
	virtual number sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const;
        virtual std::vector<number> normalised_in_box(const std::vector<number> &v);

	std::vector<number> get_abs_pos(BaseField *p);
        virtual void shift_particle(BaseField *p, std::vector<number> &amount);
	virtual void setNeighborsPeriodic(int Lx, int Ly);
	virtual number get_shear_displacement(){return delta_x;};
};

#endif /* LEESEDWARDSBOX_H_ */
