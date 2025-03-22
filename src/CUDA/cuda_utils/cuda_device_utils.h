#ifndef CUDA_DEVICE_INFO_H_
#define CUDA_DEVICE_INFO_H_

int get_device_count();
void check_device_existance(int device);
cudaDeviceProp get_current_device_prop();
cudaDeviceProp get_device_prop(int device);

#endif /* CUDA_DEVICE_INFO_H_ */
