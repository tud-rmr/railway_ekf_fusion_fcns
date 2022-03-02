function [x_hat,P,nu,S] = calcDiscreteKf(u,z,M,R,z_selector,x_0,P_0,state_transition_fcn_name,measurement_fcn_name,sample_time)
% [x_hat,P,nu,S] = calcDiscreteKf(u,z,M,R,z_selector,x_0,P_0,state_transition_fcn_name,measurement_fcn_name,sample_time)
%
%   In:
%       u                           Input vector
%       z                           Measurement vector
%       M                           Input gain matrix
%       R                           Measurement covariance matrix
%       z_selector                  Selector for valid measurments
%       x_0                         Previous state estimate
%       P_0                         Previous state estimate covariance matrix
%       state_transition_fcn_name   Name of the state transition function (this function has to provide the following arguments: in --> [x,u,sample_time], out --> [x_hat,F,Q])
%       measurement_fcn_name        Name of the measurement function (this function has to provide the following arguments: in --> [x,u,sample_time], out --> [z_hat,H])
%       sample_time                 Discrete sample time in seconds
%
%   Out:
%       x_hat                       New state estimate
%       P                           New state estimate covariance matrix
%       nu                          Measurement residuum
%       S                           Measurement residuum covariance matrix
%
%   Other m-files required: none
%   Subfunctions: isPositiveDefinite
%   MAT-files required: none
%
%   See also: none

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 26-Nov-2019

%% Initialization and checks

symmetry_round_tolerance = 1e-12;

% Checks __________________________________________________________________

if (sample_time < 0)
    error('calcDiscreteKf: Wrong sample time!');
end % if

if ~isPositiveDefinite(M,symmetry_round_tolerance)
    error('calcDiscreteKf: ''M'' is not positive semidefinite');
end % if

if ~isPositiveDefinite(R,symmetry_round_tolerance)
    error('calcDiscreteKf: ''R'' is not positive semidefinite');
end % if

if ~isPositiveDefinite(P_0,symmetry_round_tolerance)
    error('calcDiscreteKf: ''P_0'' is not positive semidefinite');
end % if

%% Calculations

% Prediction ______________________________________________________________

[x_hat,F,Q] = eval([state_transition_fcn_name,'(x_0,z,u,M,sample_time);']);
P = F*P_0*F.' + Q;
P = 1/2*(P+P.');

% Correction ______________________________________________________________

if z_selector % Valid measurement available
    [z_hat,H] = eval([measurement_fcn_name,'(x_hat,z,u,sample_time);']);
            
    nu = z - z_hat;
    
    S = H*P*H.' + R;
    S = 1/2*(S+S.');

    K = P*H'*S^(-1);

    x_hat = x_hat + K*nu;
    A = eye(length(K))-K*H;
    P = A*P*A.' + K*R*K.';
    P = 1/2*(P+P.');
else % No valid measurement available
    nu = zeros(size(z));   
    S = zeros(length(z));
end % if

%% Final checks

if ~isPositiveDefinite(Q,symmetry_round_tolerance)
    error('calcDiscreteKf: ''Q'' is not positive semidefinite');
end % if

if ~isPositiveDefinite(S,symmetry_round_tolerance)
    error('calcDiscreteKf: ''S'' is not positive semidefinite');
end % if

if ~isPositiveDefinite(P,symmetry_round_tolerance)
    error('calcDiscreteKf: ''P'' is not positive semidefinite');
end % if

end % function

%% Subfunctions

function q = isPositiveDefinite(A,tolerance)
% q = isPositiveDefinite(A,tolerance)
%
%   In:
%       A           Covariance matrix
%       tolerance   Symmetry tolerance
%
%   Out:
%       q           True if A is positive definite, otherwise q is false
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calcDiscreteKf

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 26-Nov-2019

test_sym = A - (A+A.')/2;
test_det = det(A);
if any(any(abs(test_sym) > tolerance)) || (test_det < -tolerance)
    q = false;
else
    q = true;
end % if

end

