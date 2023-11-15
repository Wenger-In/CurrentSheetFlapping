clear; close all;
plot_time = '2021-01-17 12:00:00';
Rs = 6.955e5; % [km]

%% Earth position in HCI coordinate, obtained from get_Earth_position.py
pos_earth_HCI = [41.31776271, -4.80444563, 0.98379865]; % (lon, lat, distance) in (deg, deg, AU)

%% STEREO-A position in HCI coordinate, obtained from cdaweb
pos_stereo_HCI = [345, 1.8]; % (lon, lat) in (deg, deg)

%% PSP position in HCI coordinate
% import data
spc_dir = 'E:\Research\Data\PSP\Encounter 7\';
spc_file = ['psp_swp_spc_l3i_',plot_time(1:4),plot_time(6:7),plot_time(9:10),'_v02.cdf'];
spc_path = [spc_dir,spc_file];
spc_info = spdfcdfinfo(spc_path);

% obtain PSP position
spc_Epoch = spdfcdfread(spc_path,'Variables','Epoch');
sc_pos_HCI = spdfcdfread(spc_path,'Variables','sc_pos_HCI');

num_epoch = length(spc_Epoch);
sc_pos_HCI = sc_pos_HCI / Rs;
pos_psp_HCI = sc_pos_HCI(num_epoch/2,:);

%% plot relative position
figure()
scatter(10*cosd(pos_earth_HCI(1)),10*sind(pos_earth_HCI(1))); hold on
scatter(10*cosd(pos_stereo_HCI(1)),10*sind(pos_stereo_HCI(1))); hold on
plot(sc_pos_HCI(:,1),sc_pos_HCI(:,2)); hold on
scatter(pos_psp_HCI(1),pos_psp_HCI(2));
legend();
axis equal