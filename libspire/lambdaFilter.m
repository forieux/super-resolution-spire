function value = lambdaFilter(band, lambda)
%% LAMDAFILTER - Return the value of the filter of one band
%% 
%% value = lambdaFilter(band, lambda)
%% 
%% return the value of the wavelength filter at current LAMDA. Return -1 if
%% error
%% 
%% FUNCTION CALL
%% 
%% value = lambdaFilter(band, lamda)
%% 
%% INPUT PARAMETERS
%% 
%% band -- the band 1 for 250, 2 for 260 and 3 for 520
%% 
%% lambda -- a Nx1 vector of value of lambda
  
  value = -1;
  
  if band == 1
    value = 0.521*ones(size(lambda));
  elseif band == 2
    value = 0.483*ones(size(lambda));    
  elseif band == 3
    value = 0.579*ones(size(lambda));
  else
    disp('Unknown band')
  end
    
end
