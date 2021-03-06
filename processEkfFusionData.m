function ref_processed_data = processEkfFusionData(out_filename,out_variable_name,x_hat,P_x_hat,utm_zone)
% ref_processed_data = processEkfFusionData(out_filename,out_variable_name,x_hat,P_x_hat)
% 
%   Process output data of EKF fusion
%
%   In:
%       out_filename        Name of the CSV-file and .mat-file which will be generated
%       out_variable_name   Name of variable which will be written in a .mat-file
%       x_hat               State vector fusion result from simulation output as timeseries or cell array of timeseries
%       P_x_hat             State vector covariance matrix from simulation output as timeseries or cell array of timeseries
%       utm_zone            UTM-Zone of the (x,y)-Position data in x_hat (as cell array of timeseries)
% 
%   Out:
%       sim_processed_data  Processed simulation data as timetable cell-array
%
%   Other m-files required: calcErrorEllipseParameters
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: none

%   Author: Hanno Winter
%   Date: 26-Nov-2019; Last revision: 25-March-2021

%% Init

if ~iscell(x_hat)
    x_hat = {x_hat};
end % if

if ~iscell(P_x_hat)
    P_x_hat = {P_x_hat};
end % if

sessions = find(cellfun(@(cell) ~isempty(cell),x_hat));

%% Calculations

for session_i = sessions(:)'
    
    x_hat_i = x_hat{session_i,1};
    P_x_hat_i = P_x_hat{session_i,1};
    utm_zone_i = utm_zone{session_i,1};

    fprintf('Processing EKF fusion simulation data of session %02i...',session_i);

    % Prepare data for export to csv file _________________________________

    time_unix = x_hat_i.Time;
    time_utc = datetime(time_unix,'ConvertFrom','posixtime');
    time_utc.TimeZone = 'utc';
    time_utc.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSSSSSS';    

    [ekf_fusion_lat,ekf_fusion_lon] = utm2ll( ... 
                                              x_hat_i.Data(:,1), ... 
                                              x_hat_i.Data(:,2), ... 
                                              utm_zone_i.Data, ... 
                                              'wgs84' ... 
                                            );
    sigma_latitude = sqrt(squeeze(P_x_hat_i.Data(1,1,:)));
    sigma_longitude = sqrt(squeeze(P_x_hat_i.Data(2,2,:)));
    ellip_altitude = nan(size(ekf_fusion_lat));
    sigma_ellip_altitude = nan(size(ekf_fusion_lat));

    d_track = x_hat_i.Data(:,3);
    sigma_d_track = sqrt(squeeze(P_x_hat_i.Data(3,3,:)));
    d_vehicle = cumsum([0;abs(diff(d_track))]);
    sigma_d_vehicle = sigma_d_track;
    v_vehicle = x_hat_i.Data(:,4);
    sigma_v_vehicle = sqrt(squeeze(P_x_hat_i.Data(4,4,:)));

    v_north = nan(size(v_vehicle));
    sigma_v_north = nan(size(v_vehicle));
    v_east = nan(size(v_vehicle));
    sigma_v_east = nan(size(v_vehicle));
    v_down = nan(size(v_vehicle));
    sigma_v_down = nan(size(v_vehicle));
    v_ground = abs(v_vehicle);
    sigma_v_ground = sigma_v_vehicle; 
    heading = mod(x_hat_i.Data(:,6)*180/pi,360);
    heading_selector = (-1 == sign(v_vehicle));    
    heading(heading_selector) = mod(heading(heading_selector)+180,360);
    sigma_heading = sqrt(squeeze(P_x_hat_i.Data(6,6,:)))*180/pi; 

    acc_x = x_hat_i.Data(:,5);
    acc_y = nan(size(v_ground));
    acc_z = nan(size(v_ground));

    w_x = nan(size(v_ground));
    w_y = nan(size(v_ground));
    w_z = x_hat_i.Data(:,7);

    roll = nan(size(v_ground));
    sigma_roll = nan(size(v_ground));
    pitch = nan(size(v_ground));
    sigma_pitch = nan(size(v_ground));
    yaw = x_hat_i.Data(:,6)*180/pi;
    sigma_yaw = sigma_heading;
        
    error_ellipse_confidence = 0.3934; % 0.3934 corresponds to 1-sigma 2D-error-ellipse 
    [std_major,std_minor,alpha] = calcErrorEllipseParameters(P_x_hat_i.Data(1:2,1:2,:),error_ellipse_confidence);

    ref_processed_data{session_i,1} = timetable( ...
                                                 seconds(time_unix), ...
                                                 ekf_fusion_lat, ...
                                                 ekf_fusion_lon, ...                                           
                                                 ellip_altitude, .... 
                                                 d_vehicle, ... 
                                                 d_track, ... 
                                                 v_vehicle, ... 
                                                 v_ground, ...
                                                 heading, ... 
                                                 v_north, ... 
                                                 v_east, ... 
                                                 v_down, ...     
                                                 sigma_latitude, ...
                                                 sigma_longitude, ...
                                                 sigma_ellip_altitude, ... 
                                                 sigma_d_vehicle, ... 
                                                 sigma_d_track, ... 
                                                 sigma_v_vehicle, ...
                                                 sigma_v_ground, ... 
                                                 sigma_heading, ... 
                                                 sigma_v_north, ...
                                                 sigma_v_east, ...
                                                 sigma_v_down, ...
                                                 std_major(:), ...
                                                 std_minor(:), ...
                                                 alpha(:), ...
                                                 acc_x, ...
                                                 acc_y, ...
                                                 acc_z, ...
                                                 w_x, ...
                                                 w_y, ...
                                                 w_z, ... 
                                                 roll, ...
                                                 pitch, ...
                                                 yaw, ...
                                                 sigma_roll, ... 
                                                 sigma_pitch, ... 
                                                 sigma_yaw ... 
                                               );
    ref_processed_data{session_i,1}.Properties.VariableNames = { ...
                                                                 'Latitude_deg', ... 
                                                                 'Longitude_deg', ... 
                                                                 'AltitudeEllipsoid_m', ... 
                                                                 'DistanceVehicle_m', ... 
                                                                 'DistanceTrack_m', ... 
                                                                 'VelocityVehicle_ms', ... 
                                                                 'VelocityGround_ms', ... 
                                                                 'Heading_deg', ...                                                                        
                                                                 'VelocityNorth_ms', ... 
                                                                 'VelocityEast_ms', ... 
                                                                 'VelocityDown_ms', ... 
                                                                 'LatitudeSigma_m', ... 
                                                                 'LongitudeSigma_m', ... 
                                                                 'AltitudeEllipsoidSigma_m', ...
                                                                 'DistanceVehicleSigma_m', ... 
                                                                 'DistanceTrackSigma_m', ... 
                                                                 'VelocityVehicleSigma_ms', ... 
                                                                 'VelocityGroundSigma_ms', ... 
                                                                 'HeadingSigma_deg', ... 
                                                                 'VelocityNorthSigma_ms', ... 
                                                                 'VelocityEastSigma_ms', ... 
                                                                 'VelocityDownSigma_ms', ...
                                                                 'ErrorEllipseMajor_m', ...
                                                                 'ErrorEllipseMinor_m', ...
                                                                 'ErrorEllipseOrientation_deg', ...
                                                                 'AccX_mss', ...
                                                                 'AccY_mss', ...
                                                                 'AccZ_mss', ...
                                                                 'TurnRateX_degs', ...
                                                                 'TurnRateY_degs', ...
                                                                 'TurnRateZ_degs', ...
                                                                 'Roll_deg', ...
                                                                 'Pitch_deg', ...
                                                                 'Yaw_deg', ...
                                                                 'RollSigma_deg', ...
                                                                 'PitchSigma_deg', ...
                                                                 'YawSigma_deg' ...
                                                               };

    % Finish ______________________________________________________________

    fprintf('done!\n');

end % for session_i

%% Save data 

exportToFile(ref_processed_data,out_filename,out_variable_name,'ChunkSize',500);

end % function
