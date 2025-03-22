#ifndef CUDAFORCES_H_
#define CUDAFORCES_H_

#include "../Forces/ChannelWalls.h"

#include "CUDAUtils.h"

#define CUDA_TRAP_NO_FORCE -1
#define CUDA_CHANNEL_WALLS 0

/**
 * @brief CUDA version of a ConstantRateForce.
 */
struct channel_walls{
	int type;
	c_number x, y, z;
	bool dir_as_centre;
	c_number rate;
	c_number F0;
};

void init_ChannelWalls_from_CPU(channel_walls *cuda_force, ChannelWalls *cpu_force) {
	cuda_force->type = CUDA_TRAP_CONSTANT;
	cuda_force->F0 = cpu_force->_F0;
	cuda_force->dir_as_centre = cpu_force->dir_as_centre;
	cuda_force->rate = cpu_force->_rate;
	cuda_force->x = cpu_force->_direction.x;
	cuda_force->y = cpu_force->_direction.y;
	cuda_force->z = cpu_force->_direction.z;
}

/**
 * @brief Used internally by CUDA classes to provide an inheritance-like mechanism for external forces.
 */
union CUDA_trap {
	int type;
	channel_walls channelwalls;
};

#endif /* CUDAFORCES_H_ */
