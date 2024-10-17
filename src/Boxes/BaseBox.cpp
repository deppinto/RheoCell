#include "BaseBox.h"

#include "../Fields/BaseField.h"

BaseBox::BaseBox() {

}

BaseBox::~BaseBox() {

}

//int getElementX(BaseField &p, int distX){return getElementX(p.CoM, distX)}
//int getElementX(BaseField *p, int distX){return getElementX(p->CoM, distX)}
//int getElementY(BaseField &p, int distY){return getElementY(p.CoM, distY)}
//int getElementY(BaseField *p, int distY){return getElementY(p->CoM, distY)}

std::vector<number> BaseBox::min_image(const BaseField &p, const BaseField &q) {
        return min_image(p.CoM, q.CoM);
}

std::vector<number> BaseBox::min_image(const BaseField *p, const BaseField *q) {
        return min_image(p->CoM, q->CoM);
}

number BaseBox::sqr_min_image_distance(const BaseField &p, const BaseField &q) {
        return sqr_min_image_distance(p.CoM, q.CoM);
}

number BaseBox::sqr_min_image_distance(const BaseField *p, const BaseField *q) {
        return sqr_min_image_distance(p->CoM, q->CoM);
}

