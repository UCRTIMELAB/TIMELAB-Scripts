% DERIV7_3D - 5-Tap 1st and 2nd discrete derivatives in 3D
%
% This function computes 1st and 2nd derivatives of an image using the 5-tap
% coefficients given by Farid and Simoncelli.  The results are significantly
% more accurate than MATLAB's GRADIENT function on edges that are at angles
% other than vertical or horizontal. This in turn improves gradient orientation
% estimation enormously.
%
% Usage:  [gx, gy, gxx, gyy, gxy] = derivative5(im, derivative specifiers)
%
% Arguments:
%                       im - Image to compute derivatives from.
%    derivative specifiers - A comma separated list of character strings
%                            that can be any of 'x', 'y', 'xx', 'yy' or 'xy'
%                            These can be in any order, the order of the
%                            computed output arguments will match the order
%                            of the derivative specifier strings.
% Returns:
%  Function returns requested derivatives which can be:
%     gx, gy   - 1st derivative in x and y
%     gxx, gyy - 2nd derivative in x and y
%     gxy      - 1st derivative in y of 1st derivative in x
%
%  Examples:
%    Just compute 1st derivatives in x and y
%    [gx, gy] = derivative5(im, 'x', 'y');  
%                                           
%    Compute 2nd derivative in x, 1st derivative in y and 2nd derivative in y
%    [gxx, gy, gyy] = derivative5(im, 'xx', 'y', 'yy')
%
% See also: DERIVATIVE7

% Reference: Hany Farid and Eero Simoncelli "Differentiation of Discrete
% Multi-Dimensional Signals" IEEE Trans. Image Processing. 13(4): 496-508 (2004)

% Copyright (c) 2010 Peter Kovesi
% www.peterkovesi.com
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% April 2010
%
% Adapted to 3D derivatives - Joshua Morgan 2017
% HELPFILE as above, with inclusion of z, zz, xz, and yz derivatives
% All credit to above author
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

function varargout = deriv5_3d(im, varargin)

varargin = varargin(:);
varargout = cell(size(varargin));

secondDeriv = false;
for n = 1:length(varargin)
    if length(varargin{n}) > 1
        secondDeriv = true;
        break
    end
end

if ~secondDeriv
    % 5 tap 1st derivative cofficients.  These are optimal if you are just
    % seeking the 1st deriavtives
    p = [0.037659  0.249153  0.426375  0.249153  0.037659];
    d1 =[0.109604  0.276691  0.000000 -0.276691 -0.109604];
    
    filtx2 = conv2(d1,p');
    filtx3 = convn(filtx2,permute(p,[3 1 2]));
    filty3 = permute(filtx3,[2 1 3]);
    filtz3 = permute(filtx3,[3 1 2]);
    
else
    % 5-tap 2nd derivative coefficients. The associated 1st derivative
    % coefficients are not quite as optimal as the ones above but are
    % consistent with the 2nd derivative interpolator p and thus are
    % appropriate to use if you are after both 1st and 2nd derivatives.
    p  = [0.030320  0.249724  0.439911  0.249724  0.030320];
    d1 = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
    d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
    
    filtx2 = conv2(d1,p');
    filtx3 = convn(filtx2,permute(p,[3 1 2]));
    filty3 = permute(filtx3,[2 1 3]);
    filtz3 = permute(filtx3,[3 1 2]);
    
    filt2Dx2 = conv2(d2,p');
    filt2Dx3 = convn(filt2Dx2,permute(p,[3 1 2]));
    filt2Dy3 = permute(filt2Dx3,[2 1 3]);
    filt2Dz3 = permute(filt2Dx3,[3 1 2]);
end





for n = 1:length(varargin)
    if strcmpi('x', varargin{n})
        varargout{n} = convn(im,filtx3,'same');
    elseif strcmpi('y', varargin{n})
        varargout{n} = convn(im,filty3,'same');
    elseif strcmpi('z', varargin{n})
        varargout{n} = convn(im,filtz3,'same');
    elseif strcmpi('xx', varargin{n})
        varargout{n} = convn(im,filt2Dx3,'same');
    elseif strcmpi('yy', varargin{n})
        varargout{n} = convn(im,filt2Dy3,'same');
    elseif strcmpi('zz', varargin{n})
        varargout{n} = convn(im,filt2Dz3,'same');
    elseif strcmpi('xy', varargin{n}) | strcmpi('yx', varargin{n})
        gx = convn(im,filtx3,'same');
        varargout{n} = convn(gx,filty3,'same');
    elseif strcmpi('xz', varargin{n}) | strcmpi('zx', varargin{n})
        gx = convn(im,filtx3,'same');
        varargout{n} = convn(gx,filtz3,'same');
    elseif strcmpi('yz', varargin{n}) | strcmpi('zy', varargin{n})
        gy = convn(im,filty3,'same');
        varargout{n} = convn(gy,filtz3,'same');
    else
        error(sprintf('''%s'' is an unrecognized derivative option',varargin{n}));
    end
end