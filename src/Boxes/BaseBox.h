#ifndef BASEBOX_H_
#define BASEBOX_H_

#include "../defs.h"

class BaseField;

class BaseBox {
	public:
		BaseBox();
		virtual ~BaseBox();
		virtual void init(int Lx, int Ly) = 0;

		virtual void get_settings(input_file &inp) = 0;
		virtual int getElementX(int site, int distX) = 0;
		virtual int getElementY(int site, int distY) = 0;

		virtual std::vector<number> min_image(const std::vector<number> &v1, const std::vector<number> &v2) const = 0;
		virtual std::vector<number> min_image(const BaseField &p, const BaseField &q);
		virtual std::vector<number> min_image(const BaseField *p, const BaseField *q);
		virtual number sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const = 0;
		virtual number sqr_min_image_distance(const BaseField &p, const BaseField &q);
		virtual number sqr_min_image_distance(const BaseField *p, const BaseField *q);
		virtual std::vector<number> normalised_in_box(const std::vector<number> &v) = 0;

		virtual int getXsize() = 0;
		virtual int getYsize() = 0;
		virtual int getBoxsize() = 0;

		virtual std::vector<number> get_abs_pos(BaseField *p) = 0;
		virtual void shift_particle(BaseField *p, std::vector<number> &amount) = 0;

		virtual number getWalls(int k){return 0;};

};

using BoxPtr = std::shared_ptr<BaseBox>;

#endif /* BASEBOX_H_ */
