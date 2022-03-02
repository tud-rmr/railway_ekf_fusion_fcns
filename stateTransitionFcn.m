function [x_hat,F,Q] = stateTransitionFcn(x_hat_old,z,u,M,T)
% [x_hat,F,Q] = stateTransitionFcn(x_hat_old,z,u,M,T)
%
%   In:
%       x_hat_old   Previous system state estimate ([x,y,d,v,a,theta,w]')
%       z           Measurment vector ([x,y,v_ground]')
%       u           System input vector ([a,w]')
%       M           Input covariance matrix
%       T           Sample time in seconds
%
%   Out:
%       x_hat       New state estimate
%       F           Current system matrix
%       Q           Current system covariance matrix
%
%   Other m-files required: none
%   Subfunctions: ctraModel, caModel
%   MAT-files required: none
%
%   Note:
%       At the end of this file there is a script example for the symbolic 
%       calculation of the necessary jacobians.
%   
%   See also: calcDiscreteKf

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 26-Nov-2019

%%

% If the vehicle is turning use a CTRA model
if abs(u(2)) > 1e-3

    [x_hat,F,Q] = ctraModel(x_hat_old,u,M,T);  
    
% If the vehicle is barely turning use a CA model
else
    
    [x_hat,F,Q] = caModel(x_hat_old,u,M,T);
    
end % if

end

%% Subfunctions

function [x_hat,F,Q] = ctraModel(x_hat_old,u,M,T)
% [x_hat,F,Q] = ctraModel(x_hat_old,u,M,T)
%
%   In:
%       x_hat_old   Previous system state estimate ([x,y,d,v,a,theta,w]')
%       u           System input vector ([a,w]')
%       M           Input covariance matrix
%       T           Sample time in seconds
%
%   Out:
%       x_hat       New state estimate
%       F           Current system matrix
%       Q           Current system covariance matrix
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   Note:
%       At the end of this file there is a script example for the symbolic 
%       calculation of the necessary jacobians.
%   
%   See also: stateTransitionFcn

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 26-Nov-2019

% New state estimate ______________________________________________________

u1 = u(1);
u2 = u(2);
x_hat = x_hat_old;

x     = x_hat_old(1);
y     = x_hat_old(2);
d     = x_hat_old(3);
v     = x_hat_old(4);
a     = u1;
theta = x_hat_old(6);
w     = u2;

delta_x = 1/w^2*((-v*w - a*w*T)*cos(theta+w*T) + a*sin(theta+w*T) + v*w*cos(theta) - a*sin(theta));
delta_y = 1/w^2*((v*w + a*w*T)*sin(theta+w*T) + a*cos(theta+w*T)  - v*w*sin(theta) - a*cos(theta));   
x_hat(1,1) = x + delta_x;
x_hat(2,1) = y + delta_y;
x_hat(3,1) = d + v*T + 1/2*a*T^2;
x_hat(4,1) = v + a*T;
x_hat(5,1) = a;
x_hat(6,1) = theta + w*T;
x_hat(7,1) = w;

% New linearized system matrix ____________________________________________

x     = x_hat(1);
y     = x_hat(2);
d     = x_hat(3);
v     = x_hat(4);
a     = x_hat(5);
theta = x_hat(6);
w     = x_hat(7);

F = [ ... 
      1, 0, 0, -(cos(theta + T*w) - cos(theta))/w, -(sin(theta) - sin(theta + T*w) + T*w*cos(theta + T*w))/w^2,  (a*cos(theta + T*w) - a*cos(theta) + w*sin(theta + T*w)*(v + T*a) - v*w*sin(theta))/w^2, (2*a*sin(theta) - 2*a*sin(theta + T*w) - v*w*cos(theta) + v*w*cos(theta + T*w) + T*v*w^2*sin(theta + T*w) + T^2*a*w^2*sin(theta + T*w) + 2*T*a*w*cos(theta + T*w))/w^3; ... 
      0, 1, 0,  (sin(theta + T*w) - sin(theta))/w,  (cos(theta + T*w) - cos(theta) + T*w*sin(theta + T*w))/w^2, -(a*sin(theta + T*w) - a*sin(theta) - w*cos(theta + T*w)*(v + T*a) + v*w*cos(theta))/w^2, (2*a*cos(theta) - 2*a*cos(theta + T*w) + v*w*sin(theta) - v*w*sin(theta + T*w) + T*v*w^2*cos(theta + T*w) + T^2*a*w^2*cos(theta + T*w) - 2*T*a*w*sin(theta + T*w))/w^3; ... 
      0, 0, 1,                                  T,                                                       T^2/2,                                                                                        0,                                                                                                                                                                      0; ... 
      0, 0, 0,                                  1,                                                           T,                                                                                        0,                                                                                                                                                                      0; ... 
      0, 0, 0,                                  0,                                                           1,                                                                                        0,                                                                                                                                                                      0; ... 
      0, 0, 0,                                  0,                                                           0,                                                                                        1,                                                                                                                                                                      T; ... 
      0, 0, 0,                                  0,                                                           0,                                                                                        0,                                                                                                                                                                      1  ... 
    ];

% New system covariance matrix ____________________________________________

V = [ ... 
      -(sin(theta) - sin(theta + T*u2) + T*u2*cos(theta + T*u2))/u2^2, (2*u1*sin(theta) - 2*u1*sin(theta + T*u2) - u2*v*cos(theta) + u2*v*cos(theta + T*u2) + T*u2^2*v*sin(theta + T*u2) + T^2*u1*u2^2*sin(theta + T*u2) + 2*T*u1*u2*cos(theta + T*u2))/u2^3; ... 
       (cos(theta + T*u2) - cos(theta) + T*u2*sin(theta + T*u2))/u2^2, (2*u1*cos(theta) - 2*u1*cos(theta + T*u2) + u2*v*sin(theta) - u2*v*sin(theta + T*u2) - 2*T*u1*u2*sin(theta + T*u2) + T*u2^2*v*cos(theta + T*u2) + T^2*u1*u2^2*cos(theta + T*u2))/u2^3; ... 
                                                                T^2/2,                                                                                                                                                                                     0; ... 
                                                                    T,                                                                                                                                                                                     0; ... 
                                                                    1,                                                                                                                                                                                     0; ... 
                                                                    0,                                                                                                                                                                                     T; ... 
                                                                    0,                                                                                                                                                                                     1  ... 
    ];
 
Q = V*M*V';
  
end % function

function [x_hat,F,Q] = caModel(x_hat_old,u,M,T)
% [x_hat,F,Q] = caModel(x_hat_old,u,M,T)
%
%   In:
%       x_hat_old   Previous system state estimate ([x,y,d,v,a,theta,w]')
%       u           System input vector ([a,w]')
%       M           Input covariance matrix
%       T           Sample time in seconds
%
%   Out:
%       x_hat       New state estimate
%       F           Current system matrix
%       Q           Current system covariance matrix
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   Note:
%       At the end of this file there is a script example for the symbolic 
%       calculation of the necessary jacobians.
%   
%   See also: stateTransitionFcn

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 26-Nov-2019

% New state estimate ______________________________________________________

u1 = u(1);
u2 = u(2);
x_hat = x_hat_old;

x     = x_hat_old(1);
y     = x_hat_old(2);
d     = x_hat_old(3);
v     = x_hat_old(4);
a     = u1;
theta = x_hat_old(6);
w     = u2;

delta_x = 1/2*a*T^2*sin(theta) + v*T*sin(theta);
delta_y = 1/2*a*T^2*cos(theta) + v*T*cos(theta);

x_hat(1) = x + delta_x;
x_hat(2) = y + delta_y;
x_hat(3) = d + v*T + 1/2*a*T^2;
x_hat(4) = v + a*T;
x_hat(5) = a;
x_hat(6) = theta;
x_hat(7) = 0;

% New linearized system matrix ____________________________________________

x     = x_hat(1);
y     = x_hat(2);
d     = x_hat(3);
v     = x_hat(4);
a     = x_hat(5);
theta = x_hat(6);
w     = x_hat(7);

F = [ ... 
      1, 0, 0, T*sin(theta), (T^2*sin(theta))/2,  (T*cos(theta)*(2*v + T*a))/2, 0; ... 
      0, 1, 0, T*cos(theta), (T^2*cos(theta))/2, -(T*sin(theta)*(2*v + T*a))/2, 0; ... 
      0, 0, 1,            T,              T^2/2,                             0, 0; ... 
      0, 0, 0,            1,                  T,                             0, 0; ... 
      0, 0, 0,            0,                  1,                             0, 0; ... 
      0, 0, 0,            0,                  0,                             1, 0; ... 
      0, 0, 0,            0,                  0,                             0, 0  ... 
    ];
  
% New system covariance matrix ____________________________________________

V = [ ... 
      (T^2*sin(theta))/2, 0; ... 
      (T^2*cos(theta))/2, 0; ... 
                   T^2/2, 0; ... 
                       T, 0; ... 
                       1, 0; ... 
                       0, 0; ... 
                       0, 0  ... 
    ];
Q = V*M*V';
  
end % function

%% Helper script to calculate jacobians

% % clear all
% % close all
% % clc
% 
% syms x y d v a theta w T u1 u2 real
% state_vec = [x y d v a theta w]';
% input_vec = [u1 u2];
% 
% % CTRA ____________________________________________________________________
% 
% delta_x = 1/u2^2*((-v*u2 - u1*u2*T)*cos(theta+u2*T) ... 
%           + u1*sin(theta+u2*T)                   ... 
%           + v*u2*cos(theta) - u1*sin(theta));
% delta_y = 1/u2^2*((v*u2 + u1*u2*T)*sin(theta+u2*T)  ... 
%           + u1*cos(theta+u2*T)                   ... 
%           - v*u2*sin(theta) - u1*cos(theta)); 
% x_hat(1,1) = x + delta_x;
% x_hat(2,1) = y + delta_y;
% x_hat(3,1) = d + v*T + 1/2*u1*T^2;
% x_hat(4,1) = v + u1*T;
% x_hat(5,1) = u1;
% x_hat(6,1) = theta + u2*T;
% x_hat(7,1) = u2;
% 
% F_ctra = simplify(jacobian(x_hat,state_vec))
% B_ctra = simplify(jacobian(x_hat,input_vec))
% 
% % CA ______________________________________________________________________
% 
% delta_x = 1/2*u1*T^2*sin(theta) + v*T*sin(theta);
% delta_y = 1/2*u1*T^2*cos(theta) + v*T*cos(theta);
% 
% x_hat(1) = x + delta_x;
% x_hat(2) = y + delta_y;
% x_hat(3) = d + v*T + 1/2*u1*T^2;
% x_hat(4) = v + u1*T;
% x_hat(5) = u1;
% x_hat(6) = theta;
% x_hat(7) = 0;
% 
% F_ca = simplify(jacobian(x_hat,state_vec))
% B_ca = simplify(jacobian(x_hat,input_vec))

