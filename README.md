Efficient PS Wavefield Separation for Anisotropic Media in the Spatial Domain

Codes are used for test the Efficient PS Wavefield Separation for Anisotropic Media in the Spatial Domain developed by

Qingyou Zhi, College of Geophysics, Chengdu University of Technology, Chengdu, Sichuan 610059, China.

Jiachun You, College of Geophysics, Chengdu University of Technology, Chengdu, Sichuan 610059, China. 

Haipeng Hu, College of Geophysics, Chengdu University of Technology, Chengdu, Sichuan 610059, China.

Wang QinGeng, BGP, CNPC.

Xinyu Fang, State Key Laboratory of Petroleum Resources and Engineering, China University of Petroleum (Beijing), Beijing, 102249, China.

Gang Yao, State Key Laboratory of Petroleum Resources and Engineering, China University of Petroleum (Beijing), Beijing, 102249, China. 

Correspondence author: Jiachun You (youjiachun2009@163.com)

The codes includes:
elastic2d_vtidecom_layersmodel.m:The main program file is an M file, named "elastic2d_vtidecom_layersmodel.m". The code runs the  layered model. The data is in the "layered model data" folder. 
When running, the code needs to be modified to change the location where the data is read in the main program.

ani_decomposition_wavefront_phase.p :The P file is a decoupled code, which is an embodiment of the P/S wave decomposition method provided in this paper.

derivate1_fd8.p and derivate2_fd8.p: Is the 8th order accuracy finite difference to calculate the first derivative andthe 8th order accuracy finite difference to calculate the second derivative.

generate_wavenumber.p:The file "generate_wavenumber.p" is designed to generate the wave numbers (kx, kz) in a two-dimensional grid.
meal2d.p :Description: Modified 2D effective absorbing boundary condition.

modpad2d.p: Description: pading 2D model parameter for absorbing boundary condition.

perclip.p:The purpose of perclip.p is to calculate the maximum and minimum values within the data set based on the given data and the percentage of data clipping.

read_matrix.p:The function named "read_matrix" is designed to read a floating-point matrix with n x n rows and nz columns from a file.

ZhiqingyouHessx folder is used to store the data of the running result, which corresponds to the save function in the main program
