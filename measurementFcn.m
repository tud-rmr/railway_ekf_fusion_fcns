function [z_hat,H] = measurementFcn(x_hat,z,u,T)
% [z_hat,H] = measurementFcn(x_hat,z,u,sample_time)
%
%   In:
%       x_hat           Current system state estimate
%       z               Measurment vecotr
%       u               System input vector
%       T               Sample time in seconds
%
%   Out:
%       z_hat           Current measurement estimate
%       H               Current measurement matrix
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calcDiscreteKf

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 26-Nov-2019

%%

% Positive vehicle velocity
if (x_hat(4) >= 0) 
    
    H = [  ... 
           1, 0, 0, 0, 0, 0, 0; ... 
           0, 1, 0, 0, 0, 0, 0; ... 
           0, 0, 0, 1, 0, 0, 0  ...
        ];
    
% Negative vehicle velocity
else 
    
    H = [  ... 
           1, 0, 0,  0, 0, 0, 0; ... 
           0, 1, 0,  0, 0, 0, 0; ... 
           0, 0, 0, -1, 0, 0, 0  ...
        ];
    
end % if

z_hat = H*x_hat;

end

