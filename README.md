# planck_uncertainties
This C++/Python code is extracted from a small project aiming at computing the uncertainties of the EVS structure formation model resulting from the N-D cosmological parameters distribution given by Planck 2018 results.
The approach was to solve N consecutive path optimization problems. Those "paths" were computed by evaluating the 2-sigma contours of the 2D marginalisations of the N-D cosmological parameters probability distribution. By finding the points on those curves which maximized or minimized M_max(z) (the cost function was defined using the partial derivatives of M_max(z), giving the expected mass of the most massive cluster in the redshift bin beween z and z+dz), we can estimate the uncertainties of the EVS model.
Unfortunaly, this method overestimates the final uncertainty since it only computes the parameter values over 2D marginalizations. This method could not be applied to higher order problems since computing the higher order 2-sigma contours rapidly exceeds the available computational power.
Applying other optimization constraints or Monte-Carlo methods could solve this problem.

# Results
Applying this algorithm once on the mean Planck 2018 cosmological parameter values with the "plikHM_TTTEEE_lowl_lowE_BAO" data yields the following results :
![image](https://user-images.githubusercontent.com/54234406/154732203-0e9c797b-3d8e-40ab-8765-1c9dcd956cac.png)
![image](https://user-images.githubusercontent.com/54234406/154732231-157dd42c-7df2-4374-95a7-7bd387ab628b.png)
![image](https://user-images.githubusercontent.com/54234406/154732234-59cecdb3-0922-4fb6-978d-ddee51c474d2.png)
![image](https://user-images.githubusercontent.com/54234406/154732241-04771446-3691-4c34-8dc5-64d2cb0d14e3.png)
![image](https://user-images.githubusercontent.com/54234406/154732252-56219a97-270d-47d1-86f6-f9b7ff3c2ca7.png)
