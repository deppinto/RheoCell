#ifndef SRC_BOXES_LEESEDWARDSCUBICBOX_H_
#define SRC_BOXES_LEESEDWARDSCUBICBOX_H_

#include "SquareBox.h"

class LeesEdwardsSquareBox: public SquareBox {
protected:
	number factor;

public:
	LeesEdwardsSquareBox();
	virtual ~LeesEdwardsSquareBox();

	virtual void get_settings(input_file &inp);
	virtual void init(number Lx, number Ly);

	virtual int getElementX(int site, int distX);
        virtual int getElementY(int site, int distY);

	virtual std::vector<number> min_image(const std::vector<number> &v1, const std::vector<number> &v2) const;
	virtual number sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const;
};

#endif /* SRC_BOXES_LEESEDWARDSSQUAREBOX_H_ */
