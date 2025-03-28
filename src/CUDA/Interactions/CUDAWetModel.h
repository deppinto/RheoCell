#ifndef CUDAWETMODEL_H_
#define CUDAWETMODEL_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/WetModel.h"

/**
 * @brief CUDA implementation of the {@link WetModel}.
 */

class CUDAWetModel: public CUDABaseInteraction, public WetModel {
public:
	CUDAWetModel();
	virtual ~CUDAWetModel();

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);
};

#endif /* CUDAWETMODEL_H_ */
