function [alphaNew, betaNew,  pointingNew, index, coefs, Nalpha, ...
          Nbeta] = computeIndex(alpha, beta, pointing, alpha_step, ...
                                beta_step, Nbolo, Nscan)
  %% COMPUTEINDEX - Compute index and coefs from pointing and grid step
  %%
  %% This is the computation of the extraction matrix P and the
  %% coefficient P^t P
  %%
  %% INPUTS PARAMETERS
  %%
  %% alpha, beta -- the vectors of alpha and beta to be rounded
  %%
  %% pointing -- a Nscan cell with each cell contains a 2x (Nbolo x
  %% Nsample) tab that contains the pointing information of each sample.
  %% First line is alpha, second is beta.
  %%
  %% alpha_step, beta_step -- the step to round all the position (ie new
  %% = round(old/step)*step )
    
    %% Compute the approximation simply by round
    
    alphaNew = unique(round(alpha./alpha_step)*alpha_step);
    betaNew = unique(round(beta./beta_step)*beta_step);

    Nalpha = length(alphaNew);
    Nbeta = length(betaNew);
    
    pointingNew = cell(1,Nscan);
    index = cell(1,Nscan);
    coefs = cell(1,Nscan);
    
    for iscan = 1:Nscan,

        coord_scan = pointing{iscan};
        %% Round alpha pointing
        coord_scan(1,:) = round(coord_scan(1,:)./alpha_step).* alpha_step;
        %% Round beta pointing
        coord_scan(2,:) = round(coord_scan(2,:)./beta_step).* beta_step;
        %% This the new
        pointingNew{iscan} = coord_scan;
        
        %% Find index -- Now we have a vector alpha that correspond to
        %% the discretization of our space. We have the pointing. We
        %% want to make the corresonding between the value in pointing
        %% to the index in alpha that contains this value.
        index_alpha_scan = findAllIndex(alphaNew,coord_scan(1,:));
        index_beta_scan = findAllIndex(betaNew,coord_scan(2,:));
        
	%% We work on images. We have the subindex in line column. We
        %% want the index on vector (it is identical for matlab)
        index(iscan) = {sub2ind([Nalpha Nbeta], index_alpha_scan, ...
                                index_beta_scan)};
    
        %% Computation of coefs The number of time a 'pixel' is observed
        %% and fonction of index
        coefs(iscan) = {computeCoefs(index{iscan}, Nbolo, Nalpha, ...
                                     Nbeta)};
    
    end
    
end

function index = findAllIndex(vec,value)
  %% FINDALLINDEX - return the index to build a vec from another
  %%
  %% index = findAllIndex(vec,value) return index such that value =
  %% vec(index). value must be included in vec otherwise an error is
  %% returnerd.
  %%
  %% I have try to explain this function with comment. If you find a
  %% better explaination I take it because even if it is me that have
  %% done that, I have to think 15 mins to understand again. So
  %% advantages with find\,: it is vectorized -> quick.
    
    value = value(~isnan(value));    
    [uValue , None, J] = unique(value);
    %% So we can reconstruct value with value = uValue(J).
    
    %% Now we must know where in vec we found the values of uValue.  If
    %% uValue is included in vec, the intersection "inter" is equal to
    %% uValue.
    
    [inter , None, Ivec] = intersect(uValue, vec);
    try 
        if ~isempty(find((inter - uValue) ~= 0, 1))
            fprintf(2,'Value is not included in vec\n')
            return
        end
    catch exception
        fprintf(2,'Problems during computation of index : you certainly try to find a value that is not in the vector\n')
        throw(exception)
    end
    
    %% Intersect return also inter = vec(Ivec). Since inter = uValue we
    %% have uValue = vec(Ivec). We have also value = uValue(J), so
    %% Ivec(J) are the index in vec for all element in value and value =
    %% vec(Ivec(J))

    index = Ivec(J);
    
end

function coefs = computeCoefs(index,Nbolo, Nalpha, Nbeta)
%% COMPUTECOEFS - Compute the number of time a pixel as been seen
%%
%% coefs = computeCoefs(index, Nbolo, Nalpha, Nbeta)
%% 
%% return in Nalpha x Nbeta that contains the number of time the pixel as
%% been pointed in index. Nbolo is the number of bolometer.

    indexBolo = reshape(repmat(1:Nbolo, [numel(index)/Nbolo 1]), ...
                        [numel(index) 1]);

    %% The results of repmat is a 2D matrix of numel(index) column like
    %% 
    %% | 1     | 1     | ... | 1     |
    %% | 2     | 2     | ... | 2     |
    %% | .     | .     | ... | .     |
    %% | .     | .     | ... | .     |
    %% | .     | .     | ... | .     |
    %% | Nbolo | Nbolo | ... | Nbolo |
    %% 
    %% and correspond to the bolometer index
    %% 
    %% the reshape transform this in one line ie [1 1... 2 2... N N]
    %% corresponding to index convention.
    
    coefs = sparse(index, indexBolo, 1, Nalpha*Nbeta, Nbolo);
    %% coefs is a matrix with a line per bolometer and one column for one pixel
    %% and 1 when the bolometer have see this pixel. Use sparse because there
    %% is A LOT of zero and this matrix is enormous.
    
    %% We can do the sum on bolometer and it's ok :)
    coefs = sparse(reshape(sum(coefs,2), Nalpha, Nbeta));

    %% Use parse again because there is a lot of zero also here. It's an
    %% image in sky space that indicate how many a pixel as been seen
    %% during this scan.

end

