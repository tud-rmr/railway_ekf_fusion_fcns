function [x_hat_ts,P_ts,nu_ts,S_ts,utm_zone] = runEkfFusion(gnss_data,imu_data,session_selector,Ts,Tend)
% [x_hat_ts,P_ts,nu_ts,S_ts] = runEkfFusion(gnss_data_tt,imu_data_tt,session_selector,Ts,Tend)
%
%   In:
%       gnss_data       Timetable or cell arrray of timetables of processed GNSS data
%       imu_data        Timetable or cell arrray of timetables of processed IMU data
%       session_selector  Number of sessions for which the fusion is performed. If empty all available sessions are used
%       Ts              Minimum sample time in seconds. If empty the natural sample time of the data is used.
%       Tend            Filtering end time (relative to start time of the data). If empty the simulation will be performed till the end of the data.
%
%   Out:
%       x_hat_ts        State estimate as cell array of timeseries
%       P_ts            State estimate covariance matrix as cell array of timeseries
%       nu_ts           Measurement residuum as cell array of timeseries
%       S_ts            Measurement residuum covariance matrix as cell array of timeseries
%       utm_zone        UTM-Zone of the (x,y)-Position data in x_hat_ts (as cell array of timeseries)
%
%   Other m-files required: none
%   Subfunctions: ll2utm, calcDiscreteKf
%   MAT-files required: none
%
%   See also: ll2utm, calcDiscreteKf

%   Author: Hanno Winter
%   Date: 18-Feb-2020; Last revision: 18-Nov-2020

%% Settings (edit)

% Initial state ___________________________________________________________

d_0     = 0;          % in m
theta_0 = 0*pi/180;   % in rad

% Uncertainties ___________________________________________________________

% Initial state uncertainties
sigma_0_d     = 10;          % in m
sigma_0_theta = 180*pi/180;  % in rad

% Input uncertainties
sigma_acc = 0.005*9.81;  % in m/s^2
sigma_w   = 0.05*pi/180; % in rad/s

% Fallback measurement uncertainties
sigma_x = 10; % in m 
sigma_y = 10; % in m
sigma_v = 0.5; % in m/s

% GNSS selector ___________________________________________________________

max_pdop = 3;

% Models __________________________________________________________________

state_transition_fcn_name = 'stateTransitionFcn';
measurement_fcn_name = 'measurementFcn';

%% Init

if ~iscell(gnss_data)
    gnss_data = {gnss_data};
end % if

if ~iscell(imu_data)
    imu_data = {imu_data};
end % if

gnss_sessions = find(cellfun(@(cell) ~isempty(cell),gnss_data));
imu_sessions = find(cellfun(@(cell) ~isempty(cell),imu_data));
all_sessions_selector = intersect(gnss_sessions,imu_sessions);

if ~isempty(session_selector) && ~all(ismember(session_selector,all_sessions_selector))
    impossible_sessions = session_selector(~ismember(session_selector,all_sessions_selector));
    error(['runEkfFusion: Sessions [',num2str(impossible_sessions),'] not available in either the GNSS or the IMU input data!'])
end % if

if isempty(session_selector)
    session_selector = all_sessions_selector;
end % if

%% Calculations

for session_i = session_selector(:)'
    
    fprintf('\nStarting EKF fusion for session %02i:\n',session_i);
    
    % Session data __________________________________________________________
    
    gnss_data_tt = gnss_data{session_i};
    imu_data_tt = imu_data{session_i};
    time_gnss = seconds(gnss_data_tt.Time);
    time_imu = seconds(imu_data_tt.Time);
    
    % Resample data _______________________________________________________
    
    natural_Ts = min( mean(diff(time_imu)) , mean(diff(time_gnss)) );
    
    if ~isempty(Ts)
        time_gnss = seconds(gnss_data_tt.Time);
        time_imu = seconds(imu_data_tt.Time);

        gnss_data_selector = [true;false(length(time_gnss)-1,1)];
        i_old = 1;
        for i = 2:length(time_gnss)
            ts_gnss_i = time_gnss(i)-time_gnss(i_old);
            if ts_gnss_i >= Ts
                i_old = i;
                gnss_data_selector(i-1) = true;
            end % if
        end % for i

        imu_data_selector = [true;false(length(time_imu)-1,1)];
        i_old = 1;
        for i = 2:length(time_imu)
            ts_imu_i = time_imu(i)-time_imu(i_old);
            if ts_imu_i >= Ts
                i_old = i;
                imu_data_selector(i-1) = true;
            end % if
        end % for i

        gnss_data_tt = gnss_data_tt(gnss_data_selector,:);
        imu_data_tt = imu_data_tt(imu_data_selector,:);
    end % if

    % Prepare GNSS data ___________________________________________________

    valid_gnss_data_selctor = true(size(gnss_data_tt,1),1);
    if ismember('PDOP',gnss_data_tt.Properties.VariableNames) && any(~isnan(gnss_data_tt.PDOP))
        valid_gnss_data_selctor = valid_gnss_data_selctor & (gnss_data_tt.PDOP < max_pdop);
    end % if                   
    gnss_data_tt = gnss_data_tt(valid_gnss_data_selctor,:);

    % Prepare simulation time _____________________________________________

    if isempty(Tend)
        Tend_i = min(floor(seconds([ ... 
                                     gnss_data_tt.Time(end)-gnss_data_tt.Time(1), ... 
                                     imu_data_tt.Time(end)-imu_data_tt.Time(1), ... 
                                   ])));
    else
        Tend_i = Tend;
    end % if

    time_gnss = seconds(gnss_data_tt.Time);
    time_imu = seconds(imu_data_tt.Time);
    sim_time = union(time_gnss,time_imu);
    sim_time = sim_time(sim_time<=sim_time(1)+Tend_i);
    
    sim_steps = length(sim_time);
    
    % Prepare data according to simulation time ___________________________
    
    gnss_data_tt = retime(gnss_data_tt,seconds(sim_time),'fillwithmissing');
    imu_data_tt = retime(imu_data_tt,seconds(sim_time),'fillwithmissing');
    
    if ~any(~isnan(gnss_data_tt.Latitude_deg)) || ~any(~isnan(gnss_data_tt.Longitude_deg))
        error('runEkfFusion: No valid GNSS data available! Maybe you can try another PDOP limit?');
    end % if
    
    % Prepare EKF data ____________________________________________________

    % Input data
    u = [ ... 
          imu_data_tt.AccX_mss'; ... 
          imu_data_tt.TurnRateZ_degs' * pi/180; ...
        ];
    u = fillmissing(u,'previous',2);
    M = blkdiag(sigma_acc^2,sigma_w^2);
    
    % First valid GNSS measurement determines the UTM zone used
    valid_gnss_indices = find(~isnan(gnss_data_tt.Latitude_deg));
    [~,~,initial_utm_zone] = ll2utm(gnss_data_tt.Latitude_deg(valid_gnss_indices(1)),gnss_data_tt.Longitude_deg(valid_gnss_indices(1)),'wgs84');
    
    % Output data
    [utm_x,utm_y,utm_zone_session_i] = ll2utm(gnss_data_tt.Latitude_deg,gnss_data_tt.Longitude_deg,'wgs84',initial_utm_zone);
    z = [ ... 
          utm_x'; ... 
          utm_y'; ...
          gnss_data_tt.VelocityGround_ms' ... 
        ];

    R = zeros(3,3,sim_steps);
    R(1,1,:) = gnss_data_tt.LatitudeSigma_m.^2;
    R(2,2,:) = gnss_data_tt.LongitudeSigma_m.^2;
    R(3,3,:) = gnss_data_tt.VelocityGroundSigma_ms.^2;
    z_selector = ~isnan(z(1,:)') & ~isnan(z(2,:)') & ~isnan(z(3,:)');
    if any(isnan(R(1,1,:)))
        nan_selector = squeeze(isnan(R(1,1,:)));
        R(1,1,nan_selector) = sigma_x^2;     
        if any(z_selector(nan_selector))
        %if any(nan_selector(ismember(sim_time,time_gnss)))
            warning(['ekfFusion: Using fallback uncertainty of ',num2str(sigma_x),'m for x-position data!'])
        end % if
    end % if
    if any(isnan(R(2,2,:)))
        nan_selector = squeeze(isnan(R(2,2,:)));
        R(2,2,nan_selector) = sigma_y^2;
        if any(z_selector(nan_selector))
        %if any(nan_selector(ismember(sim_time,time_gnss)))
            warning(['ekfFusion: Using fallback uncertainty of ',num2str(sigma_y),'m for y-position data!'])
        end % if
    end % if
    if any(isnan(R(3,3,:)))
        nan_selector = squeeze(isnan(R(3,3,:)));
        R(3,3,nan_selector) = sigma_v^2;
        if any(z_selector(nan_selector))
        % if any(nan_selector(ismember(sim_time,time_gnss)))
            warning(['ekfFusion: Using fallback uncertainty of ',num2str(sigma_v),'m/s for velocity data!'])
        end % if
    end % if

    % Report ______________________________________________________________
    
    if isempty(Ts)
        fprintf('\tTs=%5.3fs (inherited)\n',natural_Ts);
    elseif Ts < natural_Ts
        fprintf('\tTs=%5.3fs (entered value of Ts=%fs is not possible with this data)\n',natural_Ts,Ts);
    else
        fprintf('\tTs=%5.3fs\n',Ts);
    end % if
    fprintf('\tTend=%is\n',Tend_i);

    fprintf('Runnning EKF fusion...\n');

    % Run EKF _____________________________________________________________

    t_out = [];
    x_hat = [];
    P = [];
    nu = [];
    S =[];
    ekf_utm_zone = [];

    ekf_initalized = false;
    t_i_last_ekf_call = 1;
    for t_i = 1:sim_steps
        
        % Report EKF status _______________________________________________
        
        ekfFusionStatus(t_i-1,sim_steps,5)
        
        % Init current time step __________________________________________
        
        sim_time_i = sim_time(t_i);
        utm_zone_i = utm_zone_session_i(t_i);

        sim_time_i_last_ekf_call = sim_time(t_i_last_ekf_call);
        sample_time_i = sim_time_i-sim_time_i_last_ekf_call;
        
        u_i = u(:,t_i);
        z_i = z(:,t_i);
        R_i = R(:,:,t_i);
        
        % Data checks _____________________________________________________
                
        if isnan(u_i)
            continue
        end % if
            
        if ~isnan(z_i)
            valid_z = true;
        else
            valid_z = false;
        end % if   
        
        % Initalize EKF ___________________________________________________
        
        if ~ekf_initalized && valid_z
            x_0 = [z_i(1);z_i(2);d_0;z_i(3);u_i(1);theta_0;u_i(2)];

            sigma_0_x = R_i(1,1);
            sigma_0_y = R_i(2,2);
            sigma_0_v = R_i(3,3);     
            P_0 = blkdiag(sigma_0_x^2,sigma_0_y^2,sigma_0_d^2,sigma_0_v^2,sigma_acc^2,sigma_0_theta^2,sigma_w^2);

            ekf_initalized = true;

            continue
        elseif ~ekf_initalized
            continue
        end % if
        
        % Run EKF _________________________________________________________
        
        [x_hat_tmp,P_tmp,nu_tmp,S_tmp] = calcDiscreteKf(u_i,z_i,M,R_i,valid_z,x_0,P_0,state_transition_fcn_name,measurement_fcn_name,sample_time_i);
        x_0 = x_hat_tmp;
        P_0 = P_tmp;       
        t_i_last_ekf_call = t_i;
        
        % Output __________________________________________________________

        t_out(end+1) = sim_time_i;
        %x_hat(:,end+1) = x_hat_tmp;
        x_hat(:,end+1) = [x_hat_tmp(1:4,:);u(1,t_i);x_hat_tmp(6,:);u(2,t_i)];  
        P = cat(3,P,P_tmp);
        nu(:,end+1)= nu_tmp;
        S = cat(3,S,S_tmp);
        ekf_utm_zone(:,end+1) = utm_zone_i;

    end % for t_i
    ekfFusionStatus(t_i,sim_steps,5)

    % Finish ______________________________________________________________

    x_hat_ts{session_i,1} = timeseries(x_hat',t_out);
    P_ts{session_i,1} = timeseries(P,t_out);
    nu_ts{session_i,1} = timeseries(nu',t_out);
    S_ts{session_i,1} = timeseries(S,t_out);
    utm_zone{session_i,1} = timeseries(fillmissing(ekf_utm_zone','nearest'),t_out);    

end % for session_i

end % function

%% Helper functions

function ekfFusionStatus(t_i,t_end,step_size)
% ekfFusionStatus(t_i,t_end,step_size)
%
%   In:
%       t_i           Current step
%       t_end         Number of steps needed to finish the current task
%       step_size     Step size in percent after which the current status 
%                     should be reported. Lates after 'step_size' seconds 
%                     the status will be reported anyway.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: none

%   Author: Hanno Winter
%   Date: 18-Feb-2020; Last revision: 19-Feb-2020


persistent last_call_time progress_mem_var

% Init ____________________________________________________________________

progress = t_i/t_end*100;

current_call_time = clock;

if isempty(progress_mem_var) || (progress_mem_var > progress)
    progress_mem_var = 0;
end % if

if isempty(last_call_time)
    last_call_time = current_call_time;
end % if

progress_since_last_call = (progress - progress_mem_var);

time_since_last_call = (current_call_time - last_call_time);
time_since_last_call = time_since_last_call(4)*3600 + time_since_last_call(5)*60 + time_since_last_call(6);

% Waitbar _________________________________________________________________

if progress == 0 || progress == 100 || progress_since_last_call > step_size  || time_since_last_call > step_size
    progress_mem_var = progress;
    last_call_time = current_call_time;
    fprintf('EKF fusion status: %6.2f%%\n',progress);
end % if    

end % end function

