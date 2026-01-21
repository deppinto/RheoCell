#ifndef MOVINGWALLS_H_
#define MOVINGWALLS_H_

#include "BaseBox.h"

/**
 * @brief A channel simulation box with walls.
 */

class MovingWalls: public BaseBox {
protected:
        int sidex, sidey;
        std::vector<int> sides;
        std::vector<double> walls;
        std::vector<double> laplacian_walls;
	number lambda_wall;
	number factor;
	llint last_step = -1;
	number delta_x = 0.;

public:
        MovingWalls();
        virtual ~MovingWalls();

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

	virtual number getWalls(int k);
	virtual number getLaplacianWalls(int k);
	virtual void UpdateWalls(bool wallflag);
};

#endif /* MOVINGWALLS_H_ */
