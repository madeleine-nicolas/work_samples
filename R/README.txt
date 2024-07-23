R FUNCTIONS

GW_stream-aquif.R

	Date: August 2017

	Purpose: Implementation of analytical solution to model stream/aquifer interaction in fully saturated conditions

	Input: Observed water levels in stream and in observation borehole

	Output: Simulated water levels in observation borehole, correlation coefficient for given diffusivity value (D)
	
GW_stream-aquif_optim.R

	Date: August 2017

	Purpose: Implementation of analytical solution to model stream/aquifer interaction in fully saturated conditions 
	and optimization of parameters controlling flow

	Input: Observed water levels in stream and in all observation boreholes

	Output: Simulated water levels in observation borehole, optimal parameters

GW_stream-aquif_optim_loop.R

	Date: August 2017

	Purpose: Loops GW_stream-aquif_optim.R over all boreholes

	Input: Observed water levels in stream and in observation borehole

	Output: Simulated water levels in observation borehole, optimal parameters

	
GW_mound.R
	
	Date: June 2017

	Purpose: Implementation of Hantush analytical solution (1967) to model recharge from rectangular infiltration basin 

	Input: Observed water levels in borehole and infiltration time series

	Output: Simulated water levels in observation borehole, correlation coefficient for given K and S parameters