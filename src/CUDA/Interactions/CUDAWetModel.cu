#include "CUDAWetModel.h"
#include "CUDA_WetModel.cuh"

CUDAWetModel::CUDAWetModel() {
	_edge_compatible = true;
}

CUDAWetModel::~CUDAWetModel() {

}

void CUDAWetModel::get_settings(input_file &inp) {
	WetModel::get_settings(inp);
}

void CUDAWetModel::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	WetModel::init();

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_LJ_n, &this->_n, 3 * sizeof(int)));

	COPY_ARRAY_TO_CONSTANT(MD_sqr_rcut, this->_sqr_LJ_rcut, 3);
	COPY_ARRAY_TO_CONSTANT(MD_sqr_sigma, this->_sqr_sigma, 3);
	COPY_ARRAY_TO_CONSTANT(MD_epsilon, this->_epsilon, 3);
	COPY_ARRAY_TO_CONSTANT(MD_E_cut, this->_E_cut, 3);
}

void CUDAWetModel::compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box) {

	lj_forces
		<<<this->_launch_cfg.blocks, this->_launch_cfg.threads_per_block>>>
		(d_poss, d_forces, lists->d_matrix_neighs, lists->d_number_neighs, d_box);
	CUT_CHECK_ERROR("forces_second_step lj simple_lists error");
}
