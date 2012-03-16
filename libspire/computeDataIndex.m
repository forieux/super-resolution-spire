%% COMPUTEDATAINDEX - Compute the index of valide data
  
%% that are considere valide for the sky estimation. This data
%% correspond to data aquired during scanning at constant shift. All the
%% data are used to estimate the 1/f component. But not all for the sky.
%% This index indicate which one are used.
  
%% The hypothesis is that data for ONE array and ONE scan are in the
%% format (you can also found a reshape en a 2D tab  that is exactly the
%% same thing but arrangead visually by matlab as Nsample x Nbolo)
  
%% |----------+----------+-----+----------+--
%% | bolo1_n1 | bolo1_n2 | ... | bolo2_n1 |  
%% |----------+----------+-----+----------+--
  
%% and pointing are
  
%% |----------+----------------+-----+----------------+--
%% | bolo1_n1 | bolo1_n2_alpha | ... | bolo2_n1_alpha |  
%% |----------+----------------+-----+----------------+--
%% | bolo2_n1 | bolo1_n2_beta  | ... | bolo2_n2_beta  |  
%% |----------+----------------+-----+----------------+--
  
%% Consequently the indexData that indicate, per scan, the data that are
%% valide must be in the format
  
%% |----------+----------+-----+----------+--
%% | i1_bolo1 | i2_bolo1 | ... | i1_bolo2 |  
%% |----------+----------+-----+----------+--
  
%% and extract a sub part of the data array. It's up to you to remove
%% bound of scan, glitch or want you want.
  
%% The pointing information must now correspond to this new array. In
%% other word the pointing tab must be pointing(indexData). This is
%% pointing(indexData) that you provide to the function computeIndex
%% that compute the round pointing and index of the IMAGE that must be
%% compare to the data you consider valide.
  
%% In other word, consider that indexImage is the index of the image to
%% extract to reproduce data and make data adequation. You have to do
  
%% ||data(indexData) - convoluedSky(indexImage)||^2
  
%% TADA !!
  
%% No consider c is common part for on scan on all the bolometer and is
%% replicated so it exactly in the same format as data
  
%% |----------+----------+-----+----------+--
%% | bolo1_n1 | bolo1_n2 | ... | bolo2_n1 |  
%% |----------+----------+-----+----------+--
  
%% but the signal is the same for bolo1, bolo2 etc.
  
%% Consequently for the correlated signal, with bolometer that see the
%% sky you have to do
  
%% ||data(indexData) - convoluedSky(indexImage) - c(indexdata)||^2 
  
%% And for the blind bolometer you have to to
  
%% ||data - c||^2 
  
%% easy ?
  
%% To take into account all the scan you have to put a tab for each scan
%% in a cell.
  
%% and for the three array you have to put the cell in a bigger cell 3x1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% From now I have no idea about non valide data, I have no real data, No
%% flag information etc.... I consider my simulation provided by a
%% simulator that give me only valid data.

%% Consequently, the index return all the data. And the correlated noise is
%% without break. it's mean that there is a perfect following of the scan, no
%% transition step.
    
paramsInstrument
paramsObservation

%% Actually I know (it is my simulation) that there is 250 time sample per
%% scan. So 'c' si naturally 250*N_scan_total length. 

indexCommon250 = 1:(N_time_point*N_scan_total);

%% Since I consider there is no gap a simple reshape to get per scan

indexCommon250 = reshape(indexCommon250, [], N_scan_total)';

%% So now indexCommon250 is a 2D tab with in line the scan and in column the
% index of c that match the data at a 

%% For simulation 
indexCommon360 = indexCommon250;
indexCommon520 = indexCommon250;

indexCommon = {indexCommon250 indexCommon360 indexCommon520};

%% for iscan = 1:N_scan_total
%%     indexCommon250(iscan) = {[1:numel(data250{iscan})]};
%%     indexCommon360(iscan) = {[1:numel(data360{iscan})]};
%%     indexCommon520(iscan) = {[1:numel(data520{iscan})]};
%% end
  
%% indexCommon = {indexCommon250, indexCommon360, indexCommon520};

