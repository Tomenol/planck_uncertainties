# planck_uncertainties
This C++/Python code is extracted from a small project aiming at computing the uncertainties of the EVS structure formation model resulting from the N-D cosmological parameters distribution given by Planck 2018 results.
The approach was to solve N consecutive path optimization problems. Those "paths" were computed by evaluating the 2-sigma contours of the 2D marginalisations of the N-D cosmological parameters probability distribution. By finding the points on those curves which maximized or minimized M_max(z) (the cost function was defined using the partial derivatives of M_max(z), giving the expected mass of the most massive cluster in the redshift bin beween z and z+dz), we can estimate the uncertainties of the EVS model.
Unfortunaly, this method overestimates the final uncertainty since it only computes the parameter values over 2D marginalizations. This method could not be applied to higher order problems since computing the higher order 2-sigma contours rapidly exceeds the available computational power.
Applying other optimization constraints or Monte-Carlo methods could solve this problem.

# Results
Applying this algorithm once on the mean Planck 2018 cosmological parameter values with the "plikHM_TTTEEE_lowl_lowE_BAO" data yields the following results :
![image](https://user-images.githubusercontent.com/54234406/154732869-a92f3df4-0270-4575-941a-d6baf92343f6.png)
