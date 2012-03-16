function h = sfigure(h)
% SFIGURE - Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% See also figure
    
    if nargin >= 1 
	if ishandle(h)
            set(0, 'CurrentFigure', h);
	else
            h = figure(h);
 end
    else
	h = figure;
    end
    
end
