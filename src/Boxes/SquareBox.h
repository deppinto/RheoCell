#ifndef SQUAREBOX_H_
#define SQUAREBOX_H_

#include "BaseBox.h"

/**
 * @brief A square simulation box.
 */

class SquareBox: public BaseBox {
protected:
        int side;
        std::vector<int> sides;

public:
        SquareBox();
        virtual ~SquareBox();

	virtual void get_settings(input_file &inp);
        virtual void init(int Lx, int Ly);

	virtual int getElementX(int site, int distX);
        virtual int getElementY(int site, int distY);

        virtual int getXsize();
        virtual int getYsize();
        virtual int getBoxsize();

        virtual std::vector<number> min_image(const std::vector<number> &v1, const std::vector<number> &v2) const;
        virtual number sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const;
        virtual std::vector<number> normalised_in_box(const std::vector<number> &v);

	std::vector<number> get_abs_pos(BaseField *p);
        virtual void shift_particle(BaseField *p, std::vector<number> &amount);
};

#endif /* SQUAREBOX_H_ */
