function [x,y,z,t] = CalcParkerSolution(phi_arr,theta_arr,r_arr)
% Numerical calculation of Parker Solution
%   input:   phi_arr: Carrington Latitude array in the source surface [rad.]
%          theta_arr: Carrington Longitude array in the source surface [rad.]
%              r_arr: Radial distance array to calculate Parker Solution [m]
%   output: [x,y,z]: Magnetic field lines origined from (phi_arr,theta_arr) in Cartesian coordinate system [m]
%                 t: time needed to travel along Parker Solution [s]

    T = 1e6; % [K]
    x = zeros(length(phi_arr),length(r_arr) - 1);
    y = zeros(length(phi_arr),length(r_arr) - 1);
    z = zeros(length(phi_arr),length(r_arr) - 1);
    vel = zeros(1,length(r_arr));
    for i_phi = 1 : length(phi_arr)
        Omega = GetOmega(theta_arr(i_phi));
        phi = phi_arr(i_phi);
        theta = theta_arr(i_phi);
        t = 0;
        for i_r = 1 : length(r_arr) - 1
            V = ParkerSolution(r_arr(i_r),T);
            dr = r_arr(i_r + 1) - r_arr(i_r);
            dphi = Omega * dr / V;
            phi = phi - dphi;
            dt = dr / V;
            t = t + dt;
            x(i_phi,i_r) = r_arr(i_r) * cos(theta) * cos(phi);
            y(i_phi,i_r) = r_arr(i_r) * cos(theta) * sin(phi);
            z(i_phi,i_r) = r_arr(i_r) * sin(theta);
            vel(i_r) = V;
        end
    end
end
function omega = GetOmega(carr_lat)
% Get rotation rate of the Sun
%   input:  carr_lat: Carrington latitude [rad.]
%   output:    omega: Sun's rotation rate [rad./s]

    deg2rad = pi/180; % from degree to radian
    day2sec = 24*60*60; % from day to second

    omega1 = 14.713; omega2 = -2.396; omega3 = -1.787; % [deg./day]
    omega = (omega1 + omega2 * sin(carr_lat).^2 + omega3 * sin(carr_lat).^4) * deg2rad / day2sec;
end

function v = ParkerSolution(r,T)
% Solve Parker Solution of solar wind
%   input:  r: distance from calculate point to Sun [m]
%          vc: critical velocity [m/s]
%          rc: critical distance [m];
%   output: v: velovity of Parker Solution [m/s]
    
    kB = 1.38064852e-23; % [m^2*kg*s^-2*K^-1]
    G = 6.67259e-11; % [N*m^2*kg^-2]
    Ms = 1.989e30; % [kg]
    mp = 1.6726219e-27; % [kg]
    vc = sqrt(2 * kB * T / mp); % [m/s]
    rc = G * Ms * mp / 4 / kB / T; % [m]
    
%     func = @(v) (v/vc)^2-2*log(v/vc)-4*log(r/rc)-4*(rc/r)+3;
%     options = optimoptions('fsolve','Display','none');
    syms v
    if r<rc % v<vc while r<rc
%         v = fsolve(func,[0.5*vc],options);
        v = double(vpasolve((v/vc)^2-2*log(v/vc)-4*log(r/rc)-4*(rc/r)+3,v,[0.1*vc,0.9*vc]));
    else    % v>vc while r>rc
%         v = fsolve(func,[1.5*vc],options);
        v = double(vpasolve((v/vc)^2-2*log(v/vc)-4*log(r/rc)-4*(rc/r)+3,v,[1.1*vc,9.9*vc]));
    end
end