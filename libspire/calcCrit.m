function crit = calcCrit(object, hypers, bands, Hrond250, ...
                           Hrond360, Hrond520,  coefs250, coefs360, ...
                           coefs520, regOps, objectMean, Nalpha, ...
                           Nbeta, Norder, data, index250, index360, ...
                           index520)
% CALCQUADCRIT - Compute the value of the quadratic criterion value of object
%
% val = calcQuadCrit(object, data, hypers, bands, Hrond250, Hrond360,
%  Hrond520, index250, index360, index520, regOps, objectMean, Nalpha,
%  Nbeta, Norder)
% 
% compute the criterion value at object. Data adequation and regularization
% are quadratic. This function is indep of the correlation. A good choice is
% to make the penalization on the continous function. The correlation is
% computed in Fourier space.
%
% INPUT PARAMETERS
%
% object -- the object must be an Nalpha x Nbeta x (3*Norder) tab ordered in
% 250, 360 then 520.
%
% data -- are the data in directInvariant output convention. A three cell
% one for each band. Each cell contains data for all the scan one cell for
% each scan.
% 
% hypers -- the hypersparamter value tab of (N+1) x 3 x Norder with the number
% of line is the number of regularization operator, for each band in column,
% and the third dimension for the order. The first line hypers(1,:,1) is the
% noise precision (so indep of order) for each band.
%
% bands -- a vector of 3 with 1 to compute a band ordered in 250, 360 and
% 520. Ex : [0 1 0] compute only for 360.
% 
% Hrond* -- the transfert function for each band in the directInvariant
% convention.
%
% index* -- the index of observed pixel for each band in the directInvariant
% convention.
% 
% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
% regularization operators (a diff operator, a mean operator etc...) for
% each order in Fourier space.
% 
% objectMean -- the mean of the object law.  Must be an Nalpha x Nbeta x
% (3*Norder) tab ordered in 250, 360 then 520.
%
% Nalpha, Nbeta, Norder -- the number of alpha, beta and decomposition of
% lambda (Norder == 1 imply there is only order 0);
%
% FUNCTION CALL
%
% crit = calcQuadCrit(object, data, hypers, bands, Hrond250, Hrond360,
%  Hrond520, index250, index360, index520, regOps, objectMean, Nalpha,
%  Nbeta, Norder)

    %% Init
    paramsInstrument
    paramsObservation
  
    object = real(myifft2(object));
    
    compute250 = bands(1); compute360 = bands(2); compute520 = bands(3);    
    
    %% Data adequation
    %%% Data reproduction Hx
    output = directInvariant(object, bands, Hrond250, Hrond360, Hrond520, ...
                             index250, index360, index520, Nalpha, Nbeta, ...
                             Norder);

    
    data250 = data{1}; data360 = data{2}; data520 = data{3};
    output250 = output{1}; output360 = output{2}; output520 = output{3};
    
    adeq250 = 0; adeq360 = 0; adeq520 = 0;
    
    %%% ||y - Hx||^2
    for iscan = 1:N_scan_total
        if compute250
            adeq250 = adeq250 + sum(sum((data250{iscan} - output250{iscan}).^2));
        end
        if compute360
            adeq360 = adeq360 + sum(sum((data360{iscan} - output360{iscan}).^2));
        end
        if compute520
            adeq520 = adeq520 + sum(sum((data520{iscan} - output520{iscan}).^2));
        end
    end
    
    %%% gammaB*||y - Hx||^2
    crit = hypers(1,1,1)*adeq250 + hypers(1,2,1)*adeq360 + hypers(1,3,1)* ...
           adeq520;
    
    %% Regularization gammaX*x^tDx for each band
    if compute250
        reg250 = calcRegPart(object(:,:,1:3:end), objectMean(:,:,1:3:end), ...
                             squeeze(hypers(2:end,1,:)), regOps, Norder);
        crit = crit + reg250;
    end
    
    if compute360
        reg360 = calcRegPart(object(:,:,2:3:end), objectMean(:,:,2:3:end), ...
                             squeeze(hypers(2:end,2,:)), regOps, Norder);
        crit = crit + reg360;
    end
    
    if compute520
        reg520 = calcRegPart(object(:,:,3:3:end), objectMean(:,:,3:3:end), ...
                             squeeze(hypers(2:end,3,:)), regOps, Norder);
        crit = crit + reg520;
    end
    
end

function crit = calcRegPart(object, objectMean, hypers, regOps, Norder)
% CALCREGPART - Compute the value of the quadratic regularization criterion
%
% val = calcRegPart(object, objectMean, hypers, regOps, Norder)
% 
% compute the regularization part of the criterion of object for one
% band. This function is indep of the operator. A good choice is the
% penalization on the continus function. The correlation is computed in
% Fourier space.
%
% INPUT PARAMETERS
%
% object -- the object must be an Nalpha x Nbeta x Norder
%
% objectMean -- the mean of the object law. Must be an Nalpha x Nbeta x
% Norder
%
% hypers -- the hyperparamter value tab of N x Norder with the number of line
% is the number of regularization operator and Norder column.
% 
% regOps -- a Nalpha x Nbeta x N x Norder tab that contains the N
% regularization operators (a diff operator, a mean operator etc...) for
% each order in Fourier space
% 
% Norder -- the number of decomposition of lambda (Norder == 1 imply there
% is only order 0);
%
% FUNCTION CALL
%
% val = calcRegPart(object, objectMean, hypers, regOps, Norder)
    
    crit = 0;
        
    if size(regOps,1) == size(object,1)
        %% In fourrier space !
        for iorder = 1:Norder
            %% compute Fx
            xiorder = myfft2(object(:,:,iorder) - objectMean(:,:,iorder));
            
            for ioperator = 1:size(regOps,3)
                %% gammaI (Fx^dag) Q_I (Fx)
                contrib = real(conj(xiorder).*regOps(:,:,ioperator,iorder).* ...
                               xiorder);
                crit = crit + hypers(ioperator,iorder)*sum(contrib(:)); 
            end
            
        end
    
    else
        %% In direct space !        
        for iorder = 1:Norder
            
            im = object(:,:,iorder) - objectMean(:,:,iorder);
            
            for ioperator = 1:size(regOps,3)
                reg = regOps(:,:,ioperator,iorder);
                %% Simply x Q x where c = Qx is a convolution and xc a
                %terme wise product.
                contrib = im.*conv2(object,reg,'same');
                crit = crit + hypers(ioperator, iorder)*sum(contrib(:));
            end
            
        end
        
    end
    
end
