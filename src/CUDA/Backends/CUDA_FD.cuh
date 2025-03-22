/* System constants */
__constant__ int FD_N[1];
__constant__ float FD_dt[1];

#include "../cuda_utils/CUDA_lr_common.cuh"

// this function modifies L
// Conversion of the R matrix in _get_updated_orientation above into a quaternion simplifies into the qR quaternion used in _get_updated_orientation. This is then multiplied by the original orientation quaternion to produce the updated quaternion. Quaternion multiplication is the cross_product - dot_product and is therefore non-commutative
__device__ GPU_quat _get_updated_orientation(c_number4 &L, GPU_quat &old_o) {
	c_number norm = _module(L);
	L.x /= norm;
	L.y /= norm;
	L.z /= norm;

	c_number sintheta, costheta;
	sincos(FD_dt[0] * norm, &sintheta, &costheta);
	c_number qw = (c_number) 0.5f * sqrtf(fmaxf((c_number) 0.f, (c_number) 2.f + 2.f*costheta));
	c_number winv = (c_number)1.f /qw;
	GPU_quat R = {(c_number) 0.5f*L.x*sintheta*winv, (c_number) 0.5f*L.y*sintheta*winv, (c_number) 0.5f*L.z*sintheta*winv, qw};

	return quat_multiply(old_o, R);
}

__global__ void first_step(c_number4 *poss, GPU_quat *orientations, c_number4 *list_poss, c_number4 *vels, c_number4 *Ls, c_number4 *forces, c_number4 *torques, bool *are_lists_old) {
	if(IND >= FD_N[0]) return;

	const c_number4 F = forces[IND];

	c_number4 r = poss[IND];
	c_number4 v = vels[IND];

	v.x += F.x * (FD_dt[0] * (c_number) 0.5f);
	v.y += F.y * (FD_dt[0] * (c_number) 0.5f);
	v.z += F.z * (FD_dt[0] * (c_number) 0.5f);

	r.x += v.x * FD_dt[0];
	r.y += v.y * FD_dt[0];
	r.z += v.z * FD_dt[0];

	vels[IND] = v;
	poss[IND] = r;

	const c_number4 T = torques[IND];
	c_number4 L = Ls[IND];

	L.x += T.x * (FD_dt[0] * (c_number) 0.5f);
	L.y += T.y * (FD_dt[0] * (c_number) 0.5f);
	L.z += T.z * (FD_dt[0] * (c_number) 0.5f);

	Ls[IND] = L;


	GPU_quat qold_o = orientations[IND];
	orientations[IND] = _get_updated_orientation(L, qold_o);

	// do verlet lists need to be updated?
	if(quad_distance(r, list_poss[IND]) > FD_sqr_verlet_skin[0]) are_lists_old[0] = true;
}

__global__ void compute_molecular_coms(c_number4 *mol_coms, int *particles_to_mols, int *mol_sizes, c_number4 *poss) {
	if(IND >= FD_N[0]) {
		return;
	}

	int mol_id = particles_to_mols[IND];
	c_number4 p_contrib = poss[IND] / (c_number) mol_sizes[mol_id];

	LR_atomicAddXYZ(&(mol_coms[mol_id]), p_contrib);
}

__global__ void rescale_molecular_positions(c_number4 *mol_coms, int *particles_to_mols, c_number4 *poss, c_number4 shift_factor) {
	if(IND >= FD_N[0]) {
		return;
	}

	c_number4 ppos = poss[IND];
	c_number4 mol_com = mol_coms[particles_to_mols[IND]];
	ppos.x += mol_com.x * shift_factor.x;
	ppos.y += mol_com.y * shift_factor.y;
	ppos.z += mol_com.z * shift_factor.z;
	poss[IND] = ppos;
}

__global__ void rescale_positions(c_number4 *poss, c_number4 ratio) {
	if(IND >= FD_N[0]) {
		return;
	}
	c_number4 ppos = poss[IND];
	ppos.x *= ratio.x;
	ppos.y *= ratio.y;
	ppos.z *= ratio.z;
	poss[IND] = ppos;
}

__global__ void set_external_forces(c_number4 *poss, GPU_quat *orientations, CUDA_trap *ext_forces, c_number4 *forces, c_number4 *torques, llint step, int max_ext_forces, CUDABox *box) {
	if(IND >= FD_N[0]) return;
	// if there are no external forces then just put the force to 0
	if(max_ext_forces == 0) {
		forces[IND] = make_c_number4(0, 0, 0, 0);
		torques[IND] = make_c_number4(0, 0, 0, 0);
		return;
	}

	c_number4 ppos = poss[IND];
	c_number4 F = make_c_number4(0, 0, 0, 0);
	c_number4 T = make_c_number4(0, 0, 0, 0);

	for(int i = 0; i < max_ext_forces; i++) {
		// coalesced
		CUDA_trap extF = ext_forces[FD_N[0] * i + IND];
		switch (extF.type) {
			case CUDA_TRAP_CHANNELWALLS: {
				c_number intensity = extF.constant.F0 + step*extF.constant.rate;
				float3 dir = make_float3(extF.constant.x, extF.constant.y, extF.constant.z);
				if(extF.constant.dir_as_centre) {
					dir.x -= ppos.x;
					dir.y -= ppos.y;
					dir.z -= ppos.z;
					c_number dir_mod = _module(dir);
					dir.x /= dir_mod;
					dir.y /= dir_mod;
					dir.z /= dir_mod;
				}
				F.x += dir.x*intensity;
				F.y += dir.y*intensity;
				F.z += dir.z*intensity;
				break;
			}
			default: {
				break;
			}
		}
	}

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void second_step(c_number4 *vels, c_number4 *Ls, c_number4 *forces, c_number4 *torques) {
	if(IND >= FD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 v = vels[IND];

	v.x += F.x * FD_dt[0] * (c_number) 0.5f;
	v.y += F.y * FD_dt[0] * (c_number) 0.5f;
	v.z += F.z * FD_dt[0] * (c_number) 0.5f;
	v.w = (v.x*v.x + v.y*v.y + v.z*v.z) * (c_number) 0.5f;

	vels[IND] = v;

	c_number4 T = torques[IND];
	c_number4 L = Ls[IND];

	L.x += T.x * FD_dt[0] * (c_number) 0.5f;
	L.y += T.y * FD_dt[0] * (c_number) 0.5f;
	L.z += T.z * FD_dt[0] * (c_number) 0.5f;
	L.w = (L.x*L.x + L.y*L.y + L.z*L.z) * (c_number) 0.5f;

	Ls[IND] = L;
}
