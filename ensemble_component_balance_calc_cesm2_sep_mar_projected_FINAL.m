% CESM2-LE Ensemble analysis: 1) Calculating Arctic-average drift speeds
% for CESM2-LE from 1950-2100 (smoothed and unsmoothed options).
% 2) Calculating geostrophic (sea surface height (SSH)) and ageostrophic
% (internal stress, atmospheric stress, and oceanic stress) contributions to
% ice drift (and how the model represents the sum of these quantities to
% achieve total drift velocity). 
% 3) Calculating ageostrophic components separately (relative contribution of each). 

% Define filename components, directory paths, etc.
num_ens = 10; %Number of ensemble members being analyzed.
smooth_flag = 0; %for drift speed calculation, this flag determines if data will be smoothed. 0==unsmoothed, 1==smoothed.

% NetCDF file information: file end strings for daily (e.g.,
% 19500101-20141231*.nc) and monthly (e.g., 195001-201412*.nc) data.
hist_file_end_sep = '.19500101-20141231_arctic_sep_avg.nc'; %September
hist_file_end_mo_sep = '.195001-201412_arctic_sep_avg.nc';
fut_file_end_sep = '.20150101-21001231_arctic_sep_avg.nc';
fut_file_end_mo_sep = '.201501-210012_arctic_sep_avg.nc';
hist_file_end_mar = '.19500101-20141231_arctic_mar_avg.nc'; %March
hist_file_end_mo_mar = '.195001-201412_arctic_mar_avg.nc';
fut_file_end_mar = '.20150101-21001231_arctic_mar_avg.nc';
fut_file_end_mo_mar = '.201501-210012_arctic_mar_avg.nc';

root_dir='//wsl$/ubuntu/home/jamiewa/cesm2_lens/'; % Root directory where all data is held.
% Sub-directory by variable
in_siu_dir = 'siu_d/'; %zonal drift velocity (daily)
in_siv_dir = 'siv_d/'; %meridional drift velocity (daily)
in_salt_dir = 'salinity/'; % 4D salinity (dimensions=lat, lon, time, and level)
in_strocnx_dir = 'strocnx_d/'; %zonal ocean stress force
in_strocny_dir = 'strocny_d/'; %meridional ocean stress force
in_strairx_dir = 'strairx_d/'; %zonal atmospheric stress force
in_strairy_dir = 'strairy_d/'; %meridional atmospheric stress force
in_strintx_dir = 'strintx_d/'; %zonal internal stress
in_strinty_dir = 'strinty_d/'; %meridional internal stress
in_sithick_dir = 'sithick_d/'; %sea ice thickness
in_aice_dir = 'aice_d/'; %aggregate ice area.
in_ssh_dir = 'ssh_2/'; %ssh_mon/ is in root_dir, too.
in_PSL_dir = 'PSL/'; %sea level pressure
in_TS_dir = 'TS/'; %surface temperature
in_sst_dir = 'sst/'; %sea surface temperature
% Ensemble directory information for CESM2-LE. Each dir_ens value
% represents a different model ensemble simulation.
dir_ens{1,1} = '1011.001/';
dir_ens{2,1} = '1031.002/';
dir_ens{3,1} = '1051.003/';
dir_ens{4,1} = '1071.004/';
dir_ens{5,1} = '1091.005/';
dir_ens{6,1} = '1111.006/';
dir_ens{7,1} = '1131.007/';
dir_ens{8,1} = '1151.008/';
dir_ens{9,1} = '1171.009/';
dir_ens{10,1} = '1191.010/';
% NetCDF file front and middle information
hist_file_front = 'b.e21.BHISTsmbb.f09_g17.LE2-'; %historical (1950-2014)
fut_file_front = 'b.e21.BSSP370smbb.f09_g17.LE2-'; %future (2015-2100)
atm_file_middle = '.cam.h1.'; %Atmospheric model data (daily)
ice_file_middle = '.cice.h1.'; %sea ice model data (daily)
pop_file_middle = '.pop.h.nday1.'; %ocean model data (daily)
pop_file_middle_mo = '.pop.h.'; %ocean model data (monthly)
all_years = 1950:2100; %Years analyzed
num_years = length(all_years); %Number of years in all_years

% Import data. I need dist2coast for pre-calculated rhs terms (geo, ageo, and total).
% Dist2coast: for CESM2-LE only, used to mask out (i.e., not include) data
% points that are closer than 150km to the shoreline. This will exclude
% landfast ice and only focus on free drift.
dist2coast = ncread('/home/jamiewa/cmip6_nc_files/dist2coast_cesm2_60deg.nc','dist'); %Units: km
% Angle: used to convert grid-based orthogonal components to east-west
% versus north-south lat-lon components:
angle = ncread('C:/Users/jamwa/Downloads/cesm2_lens/angle_area_ugrid_cesm2_le.nc','ANGLE'); %Units=radians
% uarea: grid area for each grid point (units==m2).
uarea = ncread('C:/Users/jamwa/Downloads/cesm2_lens/siu_d/b.e21.BHISTsmbb.f09_g17.LE2-1151.008.cice.h1.FYarea_d.19500102-19600101_arctic.nc','uarea'); %m2

% In-variable names:
in_siu_var = in_siu_dir(1:end-1); %==siu_d, daily zonal drift velocity (m/s)
in_siv_var = in_siv_dir(1:end-1); %==siv_d, daily meridional drift velocity (m/s)
in_strairx_var_filename = in_strairx_dir(1:end-1); %==strairx_d, daily zonal atmospheric stress (N/m2)
in_strairy_var_filename = in_strairy_dir(1:end-1); %==strairy_d, daily meridional atmospheric stress
in_strocnx_var_filename = in_strocnx_dir(1:end-1); %==strocnx_d, daily zonal oceanic stress
in_strocny_var_filename = in_strocny_dir(1:end-1); %==strocny_d, daily meridional oceanic stress
in_strintx_var_filename = in_strintx_dir(1:end-1); %==strintx_d, daily zonal internal stress
in_strinty_var_filename = in_strinty_dir(1:end-1); %==strinty_d, daily meridional internal stress
in_aice_var_filename = in_aice_dir(1:end-1); %==aice_d, daily aggregated ice area. (units=None)
in_PSL_var = in_PSL_dir(1:end-1); %==PSL, mean sea level pressure (Pa)
in_sst_var = 'SST'; %Units: degrees Celsius
in_salt_var = 'SALT'; %Units: grams of salt per kilogram of ocean water.
in_TS_var = in_TS_dir(1:end-1); %==surface temperature. Units==
in_ssh_var = 'SSH_2'; %==SSH_2. Units=cm (daily sea surface height).
in_speed_var = 'speed'; % sea ice drift speed (pre-calculated). Units==m/s
in_speed_squared_var = 'speed_squared';
% The following are calculated in monthly_geo_ageo_component_calculation.m
in_geo_var = 'all_geo_annual_rhs'; % geostrophic-induced drift (pre-calculated) for March.
in_geo_var_nan = 'all_geo_annual_rhs_nan'; %geostrophic-induced drift for September, which has ice-free pixels.
in_ageo_var = 'all_ageo_annual_rhs'; % ageostrophic-induced drift (pre-calculated) for March.
in_total_var = 'all_total_annual_rhs'; %total drift (calculated from ageostrophic and geostrophic data).
in_ug_ua_var = 'ug_ua';
in_vg_va_var = 'vg_va';
in_total_rhs_speed_var = 'total_rhs_speed'; %drift speed calculated from total drift.
in_total_rhs_speed_squared_var = 'total_rhs_speed_squared';
rho_ice = 917; %kg/m3 (density of ice).
omega = 0.73e-4; %For calculating coriolis
celsius_to_kelvin = 273.15; % Conversion factor: degrees Celsius to Kelvins.

% Import latitude and longitude (ULON, ULAT). This is for the displaced
% pole grid used in CESM2's sea ice and ocean components to avoid the polar
% singularity.
lat = ncread([root_dir,in_siu_dir,dir_ens{1,1},hist_file_front,dir_ens{1,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],'ULAT');
lon = ncread([root_dir,in_siu_dir,dir_ens{1,1},hist_file_front,dir_ens{1,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],'ULON');
num_lons = size(lon,1); %lat and lon are 2D arrays, with varying longitudes stored in teh rows and varying latitudes in the columns.
num_lats = size(lat,2);

% CESM2 might use fill values for lat that are greater than 90 degrees.
% These fill values are replaced with NaN's in lon, lat, uarea, and
% dist2coast.
for row_idx=1:num_lons
    for col_idx=1:num_lats
        if lat(row_idx,col_idx)>90
            lat(row_idx,col_idx)=NaN;
            lon(row_idx,col_idx)=NaN;
            uarea(row_idx,col_idx)=NaN;
            dist2coast(row_idx,col_idx)=NaN;
        end
    end
end

% Calculate f (Coriolis parameter): to be used later when individual
% ageostrophic components are analyzed.
f = 2 * omega * sind(lat); %size==same as the latitude array.

% Initialize cell arrays for drift speed. Each cell holds a different ensemble member:
speed_all_cell_sep = cell(10,1); % September
speed_all_cell_mar = cell(10,1); % March

% Import zonal drift (siu_d) for each ensemble member (1:num_ens). This is
% done separately for historical and future data, as well as for March and
% September.
for ens_index=1:num_ens
    hist_sep_temp = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],in_speed_var);
    fut_sep_temp = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_sep],in_speed_var);
    speed_all_cell_sep{ens_index,1}(:,:,:) = cat(3,hist_sep_temp,fut_sep_temp); %Concatenate imported historical and future drift speeds.
    clear hist_sep_temp
    clear fut_sep_temp
    hist_mar_temp = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_mar],in_speed_var);
    fut_mar_temp = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_mar],in_speed_var);
    speed_all_cell_mar{ens_index,1}(:,:,:) = cat(3,hist_mar_temp,fut_mar_temp); %Concatenate imported historical and fure drift speeds.
    clear hist_mar_temp
    clear fut_mar_temp
end

% Mask out drift speed grids where dist2coast(row_idx,col_idx)<150 or
% dist2coast(row_idx,col_idx)=NaN. This is achieved by cycling through each
% latitude and each longitude for each ensemble member.
for ens_index=1:num_ens
    for row_idx=1:num_lons
        for col_idx=1:num_lats
            if dist2coast(row_idx,col_idx)<150 || isnan(dist2coast(row_idx,col_idx))
                speed_all_cell_sep{ens_index,1}(row_idx,col_idx,:) = NaN;
                speed_all_cell_mar{ens_index,1}(row_idx,col_idx,:) = NaN;
                if ens_index==1 %Because each ensemble member has the same lat, lon, and uarea information, we only NaN out values based on dist2coast thresholds once.
                    lat(row_idx,col_idx)=NaN;
                    lon(row_idx,col_idx)=NaN;
                    uarea(row_idx,col_idx)=NaN;
                end
            end
        end
    end
end
            
% Arctic speeds (not spatially-averaged). Average of all ensemble members.
% Pre-allocate NaN arrays for ensemble averages. 
sep_speed_annavg_ens = nan(num_lons,num_lats,num_years);
mar_speed_annavg_ens = nan(num_lons,num_lats,num_years);

% Calculate ensemble average by cycling through each year.
for index=1:length(all_years)
    % Create temporary 3D arrays that'll hold 2D monthly average speed info
    % for each ensemble member. This is only for year==index
    temp_speed_sep = nan(num_lons,num_lats,num_ens);
    temp_speed_mar = nan(num_lons,num_lats,num_ens);
    % Extract "index" time slice from each cell, place it in
    % temp_*(:,:,cell_index)
    for cell_index=1:num_ens %Cycling through each ensemble member to place monthly-averaged speed information into temp_speed_MONTH
        temp_speed_sep(:,:,cell_index) = speed_all_cell_sep{cell_index,1}(:,:,index);
        temp_speed_mar(:,:,cell_index) = speed_all_cell_mar{cell_index,1}(:,:,index);
    end
    % Calculate the ensemble average by taking the mean of temp_* along the
    % time dimension (NaN's are not included in the average calculation). 
    % Place the resulting average in *_ens(:,:,index).
    sep_speed_annavg_ens(:,:,index) = mean(temp_speed_sep,3,'omitnan');
    mar_speed_annavg_ens(:,:,index) = mean(temp_speed_mar,3,'omitnan');
end

%Calculate gridded sea ice drift speed trends.

% Use "area_avg_speed" function (located at the end of the script) to calculate Arctic-wide speed timeseries:
% Input variables (in order): ensemble member monthy 2D speed, grid area,
% lonigitude count, latitude count, year values, and the flag value
% determining whether the final average time series will be smoothed (if==0, not smoothed by the function. if==1, smoothed by the function).

% September
run_001_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{1,1},uarea,num_lons,num_lats,all_years,smooth_flag); % Ensemble member 1011.001
run_002_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{2,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1031.002
run_003_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{3,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1051.003
run_004_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{4,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1071.004
run_005_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{5,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1091.005
run_006_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{6,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1111.006
run_007_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{7,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1131.007
run_008_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{8,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1151.008
run_009_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{9,1},uarea,num_lons,num_lats,all_years,smooth_flag); % 1171.009
run_010_smoothed_speed_sep = area_avg_speed(speed_all_cell_sep{10,1},uarea,num_lons,num_lats,all_years,smooth_flag); % Ensemble member 1191.010

% March
run_001_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{1,1},uarea,num_lons,num_lats,all_years,smooth_flag); 
run_002_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{2,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_003_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{3,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_004_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{4,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_005_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{5,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_006_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{6,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_007_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{7,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_008_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{8,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_009_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{9,1},uarea,num_lons,num_lats,all_years,smooth_flag);
run_010_smoothed_speed_mar = area_avg_speed(speed_all_cell_mar{10,1},uarea,num_lons,num_lats,all_years,smooth_flag);

if smooth_flag==1 %Since we chop off the first and last 10 years of smoothed data (via area_avg_speed), we have to match that with the year listing (only when smooth_flag==1, though).
    all_years = all_years(11:141); %Cut off 1950-1959 and 2091-2100
end

% Place all smoothed time series in a single double array to calculate
% ensemble average speed. September:
all_runs_smoothed_speeds_sep(1,:) = run_001_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(2,:) = run_002_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(3,:) = run_003_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(4,:) = run_004_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(5,:) = run_005_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(6,:) = run_006_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(7,:) = run_007_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(8,:) = run_008_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(9,:) = run_009_smoothed_speed_sep;
all_runs_smoothed_speeds_sep(10,:) = run_010_smoothed_speed_sep;

%March
all_runs_smoothed_speeds_mar(1,:) = run_001_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(2,:) = run_002_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(3,:) = run_003_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(4,:) = run_004_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(5,:) = run_005_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(6,:) = run_006_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(7,:) = run_007_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(8,:) = run_008_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(9,:) = run_009_smoothed_speed_mar;
all_runs_smoothed_speeds_mar(10,:) = run_010_smoothed_speed_mar;

% Calculate ensemble mean: Can be smoothed or unsmoothed.
ens_smoothed_speed_sep = mean(all_runs_smoothed_speeds_sep(:,:),1,'omitnan'); %Average over all runs for each year.
ens_smoothed_speed_mar = mean(all_runs_smoothed_speeds_mar(:,:),1,'omitnan');

% % NEW: smooth the ensemble (area_avg_speed smoothes individual ensemble members. This smoothes the entire ensemble).
% ens_smoothed_speed_sep = smoothdata(ens_smoothed_speed_sep,'movmean',21);
% %Toggled off and lines below are on: non-smoothed data.
% ens_smoothed_speed_sep = ens_smoothed_speed_sep(11:141);    
% all_runs_smoothed_speeds_sep = all_runs_smoothed_speeds_sep(:,11:141);

% ens_smoothed_speed_mar = smoothdata(ens_smoothed_speed_mar,'movmean',21);
% ens_smoothed_speed_mar = ens_smoothed_speed_mar(11:141);
% all_runs_smoothed_speeds_mar = all_runs_smoothed_speeds_mar(:,11:141);
% % End NEW

% Convert drift speeds from m/s to km/day to make the plots more readable.
all_runs_smoothed_speeds_sep = all_runs_smoothed_speeds_sep * (1/1000) * (3600 * 24); %m/s * (3600 seconds per hour * 24 hours per day==days) * (1 km/1000m)
ens_smoothed_speed_sep = ens_smoothed_speed_sep * (1/1000) * (3600 * 24);
all_runs_smoothed_speeds_mar = all_runs_smoothed_speeds_mar * (1/1000) * (3600 * 24);
ens_smoothed_speed_mar = ens_smoothed_speed_mar * (1/1000) * (3600 * 24);

% Plot ensemble and individual curves. September (figure 1):
plot_001 = plot(all_years,all_runs_smoothed_speeds_sep(1,:),'color',[0.25 0.25 0.25],'LineWidth',0.5); %All ensemble member colors are grey.
% plot_001.Color = [plot_001.Color 0.5];
hold on
plot_002 = plot(all_years,all_runs_smoothed_speeds_sep(2,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_002.Color = [plot_002.Color 0.5];
plot_003 = plot(all_years,all_runs_smoothed_speeds_sep(3,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_003.Color = [plot_003.Color 0.5];
plot_004 = plot(all_years,all_runs_smoothed_speeds_sep(4,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_004.Color = [plot_004.Color 0.5];
plot_005 = plot(all_years,all_runs_smoothed_speeds_sep(5,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_005.Color = [plot_005.Color 0.5];
plot_006 = plot(all_years,all_runs_smoothed_speeds_sep(6,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_006.Color = [plot_006.Color 0.5];
plot_007 = plot(all_years,all_runs_smoothed_speeds_sep(7,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_007.Color = [plot_007.Color 0.5];
plot_008 = plot(all_years,all_runs_smoothed_speeds_sep(8,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_008.Color = [plot_008.Color 0.5];
plot_009 = plot(all_years,all_runs_smoothed_speeds_sep(9,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_009.Color = [plot_009.Color 0.5];
plot_010 = plot(all_years,all_runs_smoothed_speeds_sep(10,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_010.Color = [plot_010.Color 0.5];
plot_ens = plot(all_years,ens_smoothed_speed_sep,'k','LineWidth',4); %Ensemble average is black.
set(gca,'FontSize',24)
grid on
xlim([1945 2105])
ylim([0 15])
xlabel('Year')
ylabel('Drift Speed (km/day)')
title('Arctic-Wide September Average Sea Ice Drift Speed','FontSize',30)
hold off

% Plot ensemble and individual curves. March:
figure (2)
plot_001_mar = plot(all_years,all_runs_smoothed_speeds_mar(1,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
plot_001_mar.Color = [plot_001_mar.Color 0.5];
hold on
plot_002_mar = plot(all_years,all_runs_smoothed_speeds_mar(2,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_002_mar.Color = [plot_002_mar.Color 0.5];
plot_003_mar = plot(all_years,all_runs_smoothed_speeds_mar(3,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_003_mar.Color = [plot_003_mar.Color 0.5];
plot_004_mar = plot(all_years,all_runs_smoothed_speeds_mar(4,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_004_mar.Color = [plot_004_mar.Color 0.5];
plot_005_mar = plot(all_years,all_runs_smoothed_speeds_mar(5,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_005_mar.Color = [plot_005_mar.Color 0.5];
plot_006_mar = plot(all_years,all_runs_smoothed_speeds_mar(6,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_006_mar.Color = [plot_006_mar.Color 0.5];
plot_007_mar = plot(all_years,all_runs_smoothed_speeds_mar(7,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_007_mar.Color = [plot_007_mar.Color 0.5];
plot_008_mar = plot(all_years,all_runs_smoothed_speeds_mar(8,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_008_mar.Color = [plot_008_mar.Color 0.5];
plot_009_mar = plot(all_years,all_runs_smoothed_speeds_mar(9,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_009_mar.Color = [plot_009_mar.Color 0.5];
plot_010_mar = plot(all_years,all_runs_smoothed_speeds_mar(10,:),'color',[0.25 0.25 0.25],'LineWidth',0.5);
% plot_010_mar.Color = [plot_010_mar.Color 0.5];
plot_ens_mar = plot(all_years,ens_smoothed_speed_mar,'k','LineWidth',4);
set(gca,'FontSize',24)
grid on
xlim([1945 2105])
ylim([0 15])
xlabel('Year')
ylabel('Drift Speed (km/day)')
title('Arctic-Wide March Average Sea Ice Drift Speed','FontSize',30)
hold off

%% Dynamic steady state balance section.
% Reset all_years (they were chopped off to only include 1960-2090 in the
% previous section).
all_years = 1950:2100;
chopped_flag = 0; %If chopped_flag==1, 1960-2090. Else, 1950-2100 for entire section.

% Create cell arrays for each incoming variable. There will be a total of
% 10 cells, each cell dedicated to one ensemble member.
siu_all_cell_sep = cell(num_ens,1); % zonal drift velocity for September
siu_all_cell_mar = cell(num_ens,1); % zonal drift velocity for March
siv_all_cell_sep = cell(num_ens,1); % meridional drift velocity, September
siv_all_cell_mar = cell(num_ens,1); % meridional drift velocity, March
u_rhs_geo_all_cell_sep = cell(num_ens,1); % zonal geostrophic drift for September (calculated in monthly_geo_ageo_component_calculation.m)
u_rhs_geo_all_cell_mar = cell(num_ens,1); % zonal geostrophic drift for March
v_rhs_geo_all_cell_sep = cell(num_ens,1); % meridional geostrophic drift for September
v_rhs_geo_all_cell_mar = cell(num_ens,1); % meridional geostrophic drift for March
u_rhs_ageo_all_cell_sep = cell(num_ens,1); % zonal ageostrophic drift for September
u_rhs_ageo_all_cell_mar = cell(num_ens,1); % zonal ageostrophic drift for March
v_rhs_ageo_all_cell_sep = cell(num_ens,1); % meridional ageostrophic drift for September
v_rhs_ageo_all_cell_mar = cell(num_ens,1); % meridional ageostrophic drift for March
u_rhs_total_all_cell_sep = cell(num_ens,1); % zonal sum(geo, ageo) for September
u_rhs_total_all_cell_mar = cell(num_ens,1); % zonal sum(geo, ageo) for March
v_rhs_total_all_cell_sep = cell(num_ens,1); % meridional sum(geo, ageo) for September
v_rhs_total_all_cell_mar = cell(num_ens,1); % meridional sum(geo, ageo) for March

% Import siu for each ensemble member.
for ens_index=1:num_ens
    disp(['siu import: ens member ',dir_ens{ens_index,1}])
    % September
    siu_hist_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],in_siu_var);
    siu_fut_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_sep],in_siu_var);
    siu_all_cell_sep{ens_index,1}(:,:,:) = cat(3,siu_hist_annavg_sep,siu_fut_annavg_sep); % Concatenate historical and future drift (1950-2100)
    % March
    siu_hist_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_mar],in_siu_var);
    siu_fut_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_mar],in_siu_var);
    siu_all_cell_mar{ens_index,1}(:,:,:) = cat(3,siu_hist_annavg_mar,siu_fut_annavg_mar); % Concatenate historical and future drift (1950-2100)
end

% Import siv for each ensemble member.
for ens_index=1:num_ens
    % September
    siv_hist_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_sep],in_siv_var);
    siv_fut_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_sep],in_siv_var);
    siv_all_cell_sep{ens_index,1}(:,:,:) = cat(3,siv_hist_annavg_sep,siv_fut_annavg_sep); % Concatenate historical and future drift (1950-2100)
    % March
    siv_hist_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_mar],in_siv_var);
    siv_fut_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_mar],in_siv_var);
    siv_all_cell_mar{ens_index,1}(:,:,:) = cat(3,siv_hist_annavg_mar,siv_fut_annavg_mar); % Concatenate historical and future drift (1950-2100)
end

% Import geostrophic, ageostrophic, total rhs variables
% Zonal
% geostrophic (tilt)
for ens_index=1:num_ens
    % September
    u_rhs_geo_hist_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],in_geo_var_nan);
    u_rhs_geo_fut_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_sep],in_geo_var_nan);
    u_rhs_geo_all_cell_sep{ens_index,1}(:,:,:) = cat(3,u_rhs_geo_hist_annavg_sep,u_rhs_geo_fut_annavg_sep); % Concatenate historical and future data.
    % March
    u_rhs_geo_hist_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_mar],in_geo_var);
    u_rhs_geo_fut_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_mar],in_geo_var);
    u_rhs_geo_all_cell_mar{ens_index,1}(:,:,:) = cat(3,u_rhs_geo_hist_annavg_mar,u_rhs_geo_fut_annavg_mar); % Concatenate historical and future data.
end

% ageo (all terms added together from daily data; values calculated in monthly_geo_ageo_component_calculation.m).
for ens_index=1:num_ens
    % September
    u_rhs_ageo_hist_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],in_ageo_var);
    u_rhs_ageo_fut_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_sep],in_ageo_var);
    u_rhs_ageo_all_cell_sep{ens_index,1}(:,:,:) = cat(3,u_rhs_ageo_hist_annavg_sep,u_rhs_ageo_fut_annavg_sep); % Concatenate historical and future data.
    % March
    u_rhs_ageo_hist_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_mar],in_ageo_var);
    u_rhs_ageo_fut_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_mar],in_ageo_var);
    u_rhs_ageo_all_cell_mar{ens_index,1}(:,:,:) = cat(3,u_rhs_ageo_hist_annavg_mar,u_rhs_ageo_fut_annavg_mar); % Concatenate historical and future data.
end

% total (ageo + geo)
for ens_index=1:num_ens
    % September
    u_rhs_total_hist_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_sep],in_total_var);
    u_rhs_total_fut_annavg_sep = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_sep],in_total_var);
    u_rhs_total_all_cell_sep{ens_index,1}(:,:,:) = cat(3,u_rhs_total_hist_annavg_sep,u_rhs_total_fut_annavg_sep); % Concatenate historical and future data
    % March
    u_rhs_total_hist_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,hist_file_end_mar],in_total_var);
    u_rhs_total_fut_annavg_mar = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siu_var,fut_file_end_mar],in_total_var);
    u_rhs_total_all_cell_mar{ens_index,1}(:,:,:) = cat(3,u_rhs_total_hist_annavg_mar,u_rhs_total_fut_annavg_mar); % Concatenate historical and future data
end

% Meridional
% geostrophic (tilt)
for ens_index=1:num_ens
    % September
    v_rhs_geo_hist_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_sep],in_geo_var_nan);
    v_rhs_geo_fut_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_sep],in_geo_var_nan);
    v_rhs_geo_all_cell_sep{ens_index,1}(:,:,:) = cat(3,v_rhs_geo_hist_annavg_sep,v_rhs_geo_fut_annavg_sep); % concatenate historical and future data.
    % March
    v_rhs_geo_hist_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_mar],in_geo_var);
    v_rhs_geo_fut_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_mar],in_geo_var);
    v_rhs_geo_all_cell_mar{ens_index,1}(:,:,:) = cat(3,v_rhs_geo_hist_annavg_mar,v_rhs_geo_fut_annavg_mar); % concatenate historical and future data.
end

% ageostrophic (atm + ocn + stress)
for ens_index=1:num_ens
    % September
    v_rhs_ageo_hist_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_sep],in_ageo_var);
    v_rhs_ageo_fut_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_sep],in_ageo_var);
    v_rhs_ageo_all_cell_sep{ens_index,1}(:,:,:) = cat(3,v_rhs_ageo_hist_annavg_sep,v_rhs_ageo_fut_annavg_sep); % Concatenate historical and future data.
    % March
    v_rhs_ageo_hist_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_mar],in_ageo_var);
    v_rhs_ageo_fut_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_mar],in_ageo_var);
    v_rhs_ageo_all_cell_mar{ens_index,1}(:,:,:) = cat(3,v_rhs_ageo_hist_annavg_mar,v_rhs_ageo_fut_annavg_mar); % Concatenate historical and future data.
end

% total: geostrophic + ageostrophic
for ens_index=1:num_ens
    % September
    v_rhs_total_hist_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_sep],in_total_var);
    v_rhs_total_fut_annavg_sep = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_sep],in_total_var);
    v_rhs_total_all_cell_sep{ens_index,1}(:,:,:) = cat(3,v_rhs_total_hist_annavg_sep,v_rhs_total_fut_annavg_sep); % Concatenate historical and future data.
    % March
    v_rhs_total_hist_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,hist_file_end_mar],in_total_var);
    v_rhs_total_fut_annavg_mar = ncread([root_dir,in_siv_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_siv_var,fut_file_end_mar],in_total_var);
    v_rhs_total_all_cell_mar{ens_index,1}(:,:,:) = cat(3,v_rhs_total_hist_annavg_mar,v_rhs_total_fut_annavg_mar); % Concatenate historical and future data.
end

% Delete historical and future arrays (to save space).
clear siu_hist_annavg_sep
clear siu_fut_annavg_sep
clear siu_hist_annavg_mar
clear siu_fut_annavg_mar
clear siv_hist_annavg_sep
clear siv_fut_annavg_mar
clear u_rhs_geo_hist_annavg_sep
clear u_rhs_geo_fut_annavg_sep
clear u_rhs_geo_hist_annavg_mar
clear u_rhs_geo_fut_annavg_sep
clear v_rhs_geo_hist_annavg_sep
clear v_rhs_geo_fut_annavg_sep
clear v_rhs_geo_hist_annavg_mar
clear v_rhs_geo_fut_annavg_sep
clear u_rhs_ageo_hist_annavg_sep
clear u_rhs_ageo_fut_annavg_sep
clear u_rhs_ageo_hist_annavg_mar
clear u_rhs_ageo_fut_annavg_mar
clear v_rhs_ageo_hist_annavg_sep
clear v_rhs_ageo_fut_annavg_sep
clear v_rhs_ageo_hist_annavg_mar
clear v_rhs_ageo_fut_annavg_mar
clear u_rhs_total_hist_annavg_sep
clear u_rhs_total_fut_annavg_sep
clear u_rhs_total_hist_annavg_mar
clear u_rhs_total_fut_annavg_mar
clear v_rhs_total_hist_annavg_sep
clear v_rhs_total_fut_annavg_sep
clear v_rhs_total_hist_annavg_mar
clear v_rhs_total_fut_annavg_mar

% Now, filter out cells based on dist2coast and lat/lon designations. As
% before, variable grids corresponding to dist2coast<150 OR dist2coast==NaN are NaN'ed. 
% Finally, cycle through each ensemble member
for index=1:num_ens
    % Then, each longitude
    for row_idx=1:num_lons
        % Cycle through each latitude
        for col_idx=1:num_lats
            if dist2coast(row_idx,col_idx)<150 || isnan(dist2coast(row_idx,col_idx))
                siu_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN; 
                siu_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                if index==1 %Speed (monthly) is already ensemble-averaged, so perform dist2coast weeding only when index==1 (sep/mar_speed_annavg_ens aren't cell arrays)
                    sep_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                    mar_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                end                
            end
        end
    end
end

% Final filter: NaN out grids outside of the "Tandon domain" (see Fig. 1 in Tandon et
% la., 2018: doi:10.1029/2017JC013697)
for index=1:num_ens %Finally, cycle through each ensemble member.
    for row_idx=1:num_lons %Then, cycle through each longitude
        for col_idx=1:num_lats %First, cycle through each latitude.
            if (lon(row_idx,col_idx)>=103 && lon(row_idx,col_idx)<=236) && lat(row_idx,col_idx)<68 %Latitudes should be >=68. If not, nan out the corresponding row_idx,col_idx value in speed.
                siu_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN; 
                siu_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                if index==1 %Speed (monthly) is already ensemble-averaged, so perform dist2coast weeding only when index==1 (sep/mar_speed_annavg_ens aren't cell arrays)
                    sep_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                    mar_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                end 
            elseif (lon(row_idx,col_idx)>236 || lon(row_idx,col_idx)<103) && lat(row_idx,col_idx)<79 %Latitudes in this longitude range should be at least 79 degN
                siu_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN; 
                siu_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                if index==1 %Speed (monthly) is already ensemble-averaged, so perform dist2coast weeding only when index==1 (sep/mar_speed_annavg_ens aren't cell arrays)
                    sep_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                    mar_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                end 
            elseif isnan(lon(row_idx,col_idx)) || isnan(lat(row_idx,col_idx)) %Assuming lat/lon values are NaNs and corresponding variable values are not NaN:
                siu_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_sep{index,1}(row_idx,col_idx,:)=NaN; 
                siu_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                siv_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                u_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;           
                v_rhs_geo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_ageo_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                v_rhs_total_all_cell_mar{index,1}(row_idx,col_idx,:)=NaN;
                if index==1 %Speed (monthly) is already ensemble-averaged, so perform dist2coast weeding only when index==1 (sep/mar_speed_annavg_ens aren't cell arrays)
                    sep_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                    mar_speed_annavg_ens(row_idx,col_idx,:)=NaN;
                end 
            end
        end
    end
end

% Calculate ageostrophic values using siu/siv and our calculated
% geostrophic data (i.e., u_ageo_data=
% siu_all_cell_MON{cell_index,1}(:,:,:)-u_rhs_geo_all_cell_MON{cell_index,1}(:,:,:)).
% These will be stored in the cell arrays defined below this line in the
% code.
u_ageo_from_siu_all_cell_mar = cell(num_ens,1);
u_ageo_from_siu_all_cell_sep = cell(num_ens,1);
v_ageo_from_siv_all_cell_mar = cell(num_ens,1);
v_ageo_from_siv_all_cell_sep = cell(num_ens,1);

% Calculate the model-derived difference for ageostrophic motion (i.e., use siu_d and siv_d).
for cell_index=1:num_ens
    u_ageo_from_siu_all_cell_mar{cell_index,1}(:,:,:) = siu_all_cell_mar{cell_index,1}(:,:,:)-u_rhs_geo_all_cell_mar{cell_index,1}(:,:,:);
    u_ageo_from_siu_all_cell_sep{cell_index,1}(:,:,:) = siu_all_cell_sep{cell_index,1}(:,:,:)-u_rhs_geo_all_cell_sep{cell_index,1}(:,:,:);
    v_ageo_from_siv_all_cell_mar{cell_index,1}(:,:,:) = siv_all_cell_mar{cell_index,1}(:,:,:)-v_rhs_geo_all_cell_mar{cell_index,1}(:,:,:);
    v_ageo_from_siv_all_cell_sep{cell_index,1}(:,:,:) = siv_all_cell_sep{cell_index,1}(:,:,:)-v_rhs_geo_all_cell_sep{cell_index,1}(:,:,:);
end

% Define 3D arrays for ensemble average quantities
siu_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
siv_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
u_rhs_geo_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
u_rhs_ageo_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
u_rhs_total_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
v_rhs_geo_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
v_rhs_ageo_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
v_rhs_total_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
siu_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
siv_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
u_rhs_geo_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
u_rhs_ageo_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
u_rhs_total_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
v_rhs_geo_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
v_rhs_ageo_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
v_rhs_total_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
u_ageo_from_siu_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
u_ageo_from_siu_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
v_ageo_from_siv_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
v_ageo_from_siv_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));

% Calculate ensemble average for each year.
for index=1:length(all_years)
    % First, make temporary arrays to hold lat/lon info for each ensemble
    % member for a given year "index"
    temp_siu_sep = nan(num_lons,num_lats,num_ens);
    temp_siv_sep = nan(num_lons,num_lats,num_ens);
    temp_u_geo_sep = nan(num_lons,num_lats,num_ens);
    temp_u_ageo_sep = nan(num_lons,num_lats,num_ens);
    temp_u_total_sep = nan(num_lons,num_lats,num_ens);
    temp_v_geo_sep = nan(num_lons,num_lats,num_ens);
    temp_v_ageo_sep = nan(num_lons,num_lats,num_ens);
    temp_v_total_sep = nan(num_lons,num_lats,num_ens);
    temp_siu_mar = nan(num_lons,num_lats,num_ens);
    temp_siv_mar = nan(num_lons,num_lats,num_ens);
    temp_u_geo_mar = nan(num_lons,num_lats,num_ens);
    temp_u_ageo_mar = nan(num_lons,num_lats,num_ens);
    temp_u_total_mar = nan(num_lons,num_lats,num_ens);
    temp_v_geo_mar = nan(num_lons,num_lats,num_ens);
    temp_v_ageo_mar = nan(num_lons,num_lats,num_ens);
    temp_v_total_mar = nan(num_lons,num_lats,num_ens);
    temp_u_ageo_from_siu_mar = nan(num_lons,num_lats,num_ens);
    temp_u_ageo_from_siu_sep = nan(num_lons,num_lats,num_ens);
    temp_v_ageo_from_siv_mar = nan(num_lons,num_lats,num_ens);
    temp_v_ageo_from_siv_sep = nan(num_lons,num_lats,num_ens);
    
    % Extract "index" time slice from each VAR_all_cell_MONTH, place it in
    % temp_*(:,:,cell_index)
    for cell_index=1:num_ens
        temp_siu_sep(:,:,cell_index) = siu_all_cell_sep{cell_index,1}(:,:,index);
        temp_siv_sep(:,:,cell_index) = siv_all_cell_sep{cell_index,1}(:,:,index);
        temp_u_geo_sep(:,:,cell_index) = u_rhs_geo_all_cell_sep{cell_index,1}(:,:,index);
        temp_u_ageo_sep(:,:,cell_index) = u_rhs_ageo_all_cell_sep{cell_index,1}(:,:,index);
        temp_u_total_sep(:,:,cell_index) = u_rhs_total_all_cell_sep{cell_index,1}(:,:,index);
        temp_v_geo_sep(:,:,cell_index) = v_rhs_geo_all_cell_sep{cell_index,1}(:,:,index);
        temp_v_ageo_sep(:,:,cell_index) = v_rhs_ageo_all_cell_sep{cell_index,1}(:,:,index);
        temp_v_total_sep(:,:,cell_index) = v_rhs_total_all_cell_sep{cell_index,1}(:,:,index);
        temp_siu_mar(:,:,cell_index) = siu_all_cell_mar{cell_index,1}(:,:,index);
        temp_siv_mar(:,:,cell_index) = siv_all_cell_mar{cell_index,1}(:,:,index);
        temp_u_geo_mar(:,:,cell_index) = u_rhs_geo_all_cell_mar{cell_index,1}(:,:,index);
        temp_u_ageo_mar(:,:,cell_index) = u_rhs_ageo_all_cell_mar{cell_index,1}(:,:,index);
        temp_u_total_mar(:,:,cell_index) = u_rhs_total_all_cell_mar{cell_index,1}(:,:,index);
        temp_v_geo_mar(:,:,cell_index) = v_rhs_geo_all_cell_mar{cell_index,1}(:,:,index);
        temp_v_ageo_mar(:,:,cell_index) = v_rhs_ageo_all_cell_mar{cell_index,1}(:,:,index);
        temp_v_total_mar(:,:,cell_index) = v_rhs_total_all_cell_mar{cell_index,1}(:,:,index);
        temp_u_ageo_from_siu_mar(:,:,cell_index) = u_ageo_from_siu_all_cell_mar{cell_index,1}(:,:,index);
        temp_u_ageo_from_siu_sep(:,:,cell_index) = u_ageo_from_siu_all_cell_sep{cell_index,1}(:,:,index);
        temp_v_ageo_from_siv_mar(:,:,cell_index) = v_ageo_from_siv_all_cell_mar{cell_index,1}(:,:,index);
        temp_v_ageo_from_siv_sep(:,:,cell_index) = v_ageo_from_siv_all_cell_sep{cell_index,1}(:,:,index);
    end
    % Calculate the ensemble average by taking the mean of temp_* along the
    % time dimension. Place the resulting average in *annavg_ens_MONTH(:,:,index).
    siu_all_annavg_ens_sep(:,:,index) = mean(temp_siu_sep,3,'omitnan');
    siv_all_annavg_ens_sep(:,:,index) = mean(temp_siv_sep,3,'omitnan');
    u_rhs_geo_all_annavg_ens_sep(:,:,index) = mean(temp_u_geo_sep,3,'omitnan');
    u_rhs_ageo_all_annavg_ens_sep(:,:,index) = mean(temp_u_ageo_sep,3,'omitnan');
    u_rhs_total_all_annavg_ens_sep(:,:,index) = mean(temp_u_total_sep,3,'omitnan');
    v_rhs_geo_all_annavg_ens_sep(:,:,index) = mean(temp_v_geo_sep,3,'omitnan');
    v_rhs_ageo_all_annavg_ens_sep(:,:,index) = mean(temp_v_ageo_sep,3,'omitnan');
    v_rhs_total_all_annavg_ens_sep(:,:,index) = mean(temp_v_total_sep,3,'omitnan');
    siu_all_annavg_ens_mar(:,:,index) = mean(temp_siu_mar,3,'omitnan');
    siv_all_annavg_ens_mar(:,:,index) = mean(temp_siv_mar,3,'omitnan');
    u_rhs_geo_all_annavg_ens_mar(:,:,index) = mean(temp_u_geo_mar,3,'omitnan');
    u_rhs_ageo_all_annavg_ens_mar(:,:,index) = mean(temp_u_ageo_mar,3,'omitnan');
    u_rhs_total_all_annavg_ens_mar(:,:,index) = mean(temp_u_total_mar,3,'omitnan');
    v_rhs_geo_all_annavg_ens_mar(:,:,index) = mean(temp_v_geo_mar,3,'omitnan');
    v_rhs_ageo_all_annavg_ens_mar(:,:,index) = mean(temp_v_ageo_mar,3,'omitnan');
    v_rhs_total_all_annavg_ens_mar(:,:,index) = mean(temp_v_total_mar,3,'omitnan');
    u_ageo_from_siu_all_annavg_ens_mar(:,:,index) = mean(temp_u_ageo_from_siu_mar,3,'omitnan');
    u_ageo_from_siu_all_annavg_ens_sep(:,:,index) = mean(temp_u_ageo_from_siu_sep,3,'omitnan');
    v_ageo_from_siv_all_annavg_ens_mar(:,:,index) = mean(temp_v_ageo_from_siv_mar,3,'omitnan');
    v_ageo_from_siv_all_annavg_ens_sep(:,:,index) = mean(temp_v_ageo_from_siv_sep,3,'omitnan');
end

% Calculate All/High and Low Years climatologies (indices 1-21 for
% Steady Years, 66-85 for Low)
high_years_indices = 1:21; %1950-1970
low_years_indices = 66:86; %2015-2035

% "orig" in the variable name means that the data are on their native
% displaced pole grid.
siu_all_high_ens_clim_orig_mar = mean(siu_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan'); %March (analogue to increasing speed years)
siu_low_ens_clim_orig_sep = mean(siu_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan'); %September (analogue to decreasing speed years)
siu_low_ens_clim_orig_mar = mean(siu_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan'); %March (projected years; 2025-2090/2100, like September)
siv_all_high_ens_clim_orig_mar = mean(siv_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
siv_low_ens_clim_orig_sep = mean(siv_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
siv_low_ens_clim_orig_mar = mean(siv_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
u_rhs_geo_all_high_ens_clim_orig_mar = mean(u_rhs_geo_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
u_rhs_geo_low_ens_clim_orig_sep = mean(u_rhs_geo_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
u_rhs_geo_low_ens_clim_orig_mar = mean(u_rhs_geo_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
u_rhs_ageo_all_high_ens_clim_orig_mar = mean(u_rhs_ageo_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
u_rhs_ageo_low_ens_clim_orig_sep = mean(u_rhs_ageo_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
u_rhs_ageo_low_ens_clim_orig_mar = mean(u_rhs_ageo_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
u_rhs_total_all_high_ens_clim_orig_mar = mean(u_rhs_total_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
u_rhs_total_low_ens_clim_orig_sep = mean(u_rhs_total_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
u_rhs_total_low_ens_clim_orig_mar = mean(u_rhs_total_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
v_rhs_geo_all_high_ens_clim_orig_mar = mean(v_rhs_geo_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
v_rhs_geo_low_ens_clim_orig_sep = mean(v_rhs_geo_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
v_rhs_geo_low_ens_clim_orig_mar = mean(v_rhs_geo_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
v_rhs_ageo_all_high_ens_clim_orig_mar = mean(v_rhs_ageo_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
v_rhs_ageo_low_ens_clim_orig_sep = mean(v_rhs_ageo_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
v_rhs_ageo_low_ens_clim_orig_mar = mean(v_rhs_ageo_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
v_rhs_total_all_high_ens_clim_orig_mar = mean(v_rhs_total_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
v_rhs_total_low_ens_clim_orig_sep = mean(v_rhs_total_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
v_rhs_total_low_ens_clim_orig_mar = mean(v_rhs_total_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
u_ageo_from_siu_all_high_ens_clim_orig_mar = mean(u_ageo_from_siu_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
u_ageo_from_siu_low_ens_clim_orig_mar = mean(u_ageo_from_siu_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
u_ageo_from_siu_low_ens_clim_orig_sep = mean(u_ageo_from_siu_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
v_ageo_from_siv_all_high_ens_clim_orig_mar = mean(v_ageo_from_siv_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
v_ageo_from_siv_low_ens_clim_orig_mar = mean(v_ageo_from_siv_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
v_ageo_from_siv_low_ens_clim_orig_sep = mean(v_ageo_from_siv_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');

if chopped_flag==1 %Get rid of 1950-1959, 2091-2100
    all_years = 1960:2090;
    all_years = all_years-1949; %Convert from actual years to year indices (e.g., all_years=1960 changes to all_years=11, or the eleventh element of the original time sequence).

    % Redefine ensemble average variables by cutting off the first and last
    % 10 years.
    siu_all_annavg_ens_sep = siu_all_annavg_ens_sep(:,:,all_years);
    siv_all_annavg_ens_sep = siv_all_annavg_ens_sep(:,:,all_years);
    u_rhs_geo_all_annavg_ens_sep = u_rhs_geo_all_annavg_ens_sep(:,:,all_years);
    u_rhs_ageo_all_annavg_ens_sep = u_rhs_ageo_all_annavg_ens_sep(:,:,all_years);
    u_rhs_total_all_annavg_ens_sep = u_rhs_total_all_annavg_ens_sep(:,:,all_years);
    v_rhs_geo_all_annavg_ens_sep = v_rhs_geo_all_annavg_ens_sep(:,:,all_years);
    v_rhs_ageo_all_annavg_ens_sep = v_rhs_ageo_all_annavg_ens_sep(:,:,all_years);
    v_rhs_total_all_annavg_ens_sep = v_rhs_total_all_annavg_ens_sep(:,:,all_years);
    siu_all_annavg_ens_mar = siu_all_annavg_ens_mar(:,:,all_years);
    siv_all_annavg_ens_mar = siv_all_annavg_ens_mar(:,:,all_years);
    u_rhs_geo_all_annavg_ens_mar = u_rhs_geo_all_annavg_ens_mar(:,:,all_years);
    u_rhs_ageo_all_annavg_ens_mar = u_rhs_ageo_all_annavg_ens_mar(:,:,all_years);
    u_rhs_total_all_annavg_ens_mar = u_rhs_total_all_annavg_ens_mar(:,:,all_years);
    v_rhs_geo_all_annavg_ens_mar = v_rhs_geo_all_annavg_ens_mar(:,:,all_years);
    v_rhs_ageo_all_annavg_ens_mar = v_rhs_ageo_all_annavg_ens_mar(:,:,all_years);
    v_rhs_total_all_annavg_ens_mar = v_rhs_total_all_annavg_ens_mar(:,:,all_years);

    all_years = all_years-10; % Since the data has been chopped, we want 1960 to correspond to an index of all_years(1)=1, not all_years(1)=11. Therefore, we subtract 10 from all_years.

    % From this point foward, we'll only look at ensemble trends, etc.
    % Extract low and high years from *_ens for each variable. 
    high_years = 1960:2024; %Assuming the switch happens at 2024.
    high_years = high_years-1959; %The first index==1 now.
    low_years = 2025:2090; %OR 2044-2090, if include flattening out area.
    low_years = low_years-1959;
else %Years aren't chopped off the beginning and end of the time series.
    all_years = all_years-1949; %Convert actual year values into indices (i.e., all_years(1)=1950 becomes all_years(1)=1)
    
    high_years = 1950:2024; %March historical ONLY
    high_years = high_years-1949;
    low_years = 2025:2100; %Projected, both September and March.
    low_years = low_years-1949;
end

% Extract appropriate data for each year chunk and each ensemble variable:
% siu_d
siu_low_annavg_ens_sep = siu_all_annavg_ens_sep(:,:,low_years); %SEPTEMBER, PROJECTED YEARS
siu_low_annavg_ens_mar = siu_all_annavg_ens_mar(:,:,low_years); %MARCH, PROJECTED YEARS
siu_high_annavg_ens_mar = siu_all_annavg_ens_mar(:,:,high_years); %MARCH, HISTORICAL YEARS
% siv_d
siv_low_annavg_ens_sep = siv_all_annavg_ens_sep(:,:,low_years);
siv_low_annavg_ens_mar = siv_all_annavg_ens_mar(:,:,low_years); 
siv_high_annavg_ens_mar = siv_all_annavg_ens_mar(:,:,high_years);
% u_rhs_geo
u_rhs_geo_low_annavg_ens_sep = u_rhs_geo_all_annavg_ens_sep(:,:,low_years);
u_rhs_geo_low_annavg_ens_mar = u_rhs_geo_all_annavg_ens_mar(:,:,low_years); 
u_rhs_geo_high_annavg_ens_mar = u_rhs_geo_all_annavg_ens_mar(:,:,high_years);
% u_rhs_ageo
u_rhs_ageo_low_annavg_ens_sep = u_rhs_ageo_all_annavg_ens_sep(:,:,low_years);
u_rhs_ageo_low_annavg_ens_mar = u_rhs_ageo_all_annavg_ens_mar(:,:,low_years); 
u_rhs_ageo_high_annavg_ens_mar = u_rhs_ageo_all_annavg_ens_mar(:,:,high_years);
% u_rhs_total
u_rhs_total_low_annavg_ens_sep = u_rhs_total_all_annavg_ens_sep(:,:,low_years);
u_rhs_total_low_annavg_ens_mar = u_rhs_total_all_annavg_ens_mar(:,:,low_years); 
u_rhs_total_high_annavg_ens_mar = u_rhs_total_all_annavg_ens_mar(:,:,high_years);
% v_rhs_geo
v_rhs_geo_low_annavg_ens_sep = v_rhs_geo_all_annavg_ens_sep(:,:,low_years);
v_rhs_geo_low_annavg_ens_mar = v_rhs_geo_all_annavg_ens_mar(:,:,low_years); 
v_rhs_geo_high_annavg_ens_mar = v_rhs_geo_all_annavg_ens_mar(:,:,high_years);
% v_rhs_ageo
v_rhs_ageo_low_annavg_ens_sep = v_rhs_ageo_all_annavg_ens_sep(:,:,low_years);
v_rhs_ageo_low_annavg_ens_mar = v_rhs_ageo_all_annavg_ens_mar(:,:,low_years); 
v_rhs_ageo_high_annavg_ens_mar = v_rhs_ageo_all_annavg_ens_mar(:,:,high_years);
% v_rhs_total
v_rhs_total_low_annavg_ens_sep = v_rhs_total_all_annavg_ens_sep(:,:,low_years);
v_rhs_total_low_annavg_ens_mar = v_rhs_total_all_annavg_ens_mar(:,:,low_years); 
v_rhs_total_high_annavg_ens_mar = v_rhs_total_all_annavg_ens_mar(:,:,high_years);
% Speed (by grid)
sep_speed_low_annavg_ens = sep_speed_annavg_ens(:,:,low_years);
mar_speed_low_annavg_ens = mar_speed_annavg_ens(:,:,low_years);
% siu-geo
u_ageo_from_siu_low_annavg_ens_mar = u_ageo_from_siu_all_annavg_ens_mar(:,:,low_years);
u_ageo_from_siu_low_annavg_ens_sep = u_ageo_from_siu_all_annavg_ens_sep(:,:,low_years);
% siv-geo
v_ageo_from_siv_low_annavg_ens_mar = v_ageo_from_siv_all_annavg_ens_mar(:,:,low_years);
v_ageo_from_siv_low_annavg_ens_sep = v_ageo_from_siv_all_annavg_ens_sep(:,:,low_years);

% In September, some lat/lon combinations may only have a couple of non-ice
% free years. We delete/NaN out these grid points to avoid erroneous trend
% calculations; pixels with only 5 years of non-NaN
% entries will be completely NaN'ed before trends are calculated.

for row_idx=1:num_lons %Finally, cycle by longitude
    for col_idx=1:num_lats % Then, cycle by latitude
        num_non_nans = 0;
        for time_index=1:length(low_years) % First, cycle through each of the decreasing speed projected years.
            if ~isnan(siu_low_annavg_ens_sep(row_idx,col_idx,time_index)) %Determine whether a given lat/lon grid for year=time_index is non-NaN. 
                num_non_nans = num_non_nans + 1; %If it's non-NaN, then we will count it towards values we'll use to calculate trends.
            end
        end
        if num_non_nans<=5 && col_idx>5 && row_idx>5 % If there are only 5 or less non-NaN years, we NaN these out now (for a given lat/lon point over all years).
            for time_index=1:length(low_years) % Iterate through each year and NaN out non-NaN values
                if ~isnan(siu_low_annavg_ens_sep(row_idx,col_idx,time_index))
                    siu_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    siv_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_rhs_geo_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_rhs_geo_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_rhs_ageo_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_rhs_ageo_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_rhs_total_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_rhs_total_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_ageo_from_siu_low_annavg_ens_mar(row_idx,col_idx,time_index) = NaN;
                    u_ageo_from_siu_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_ageo_from_siv_low_annavg_ens_mar(row_idx,col_idx,time_index) = NaN;
                    v_ageo_from_siv_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                end
            end
        end
    end
end


% Calculate trends:
% Transpose high, all, and low years for ease of calculation: 
high_years = high_years';
low_years = low_years';
all_years = all_years';

% Trend calculation (angle correction will be done after this)
% trend_calc function is at the bottom of this script. Inputs (in order):
% ensemble average variable, latitude coordinates, longitude coordinates,
% applicable year chunks.
siu_all_ens_trend_mar_orig = trend_calc(siu_all_annavg_ens_mar,lat,lon,all_years); %19(5-6)0-2090 MARCH
siu_low_ens_trend_sep_orig = trend_calc(siu_low_annavg_ens_sep,lat,lon,low_years); %2025-2(0,1)(9,0)0 SEPT
siu_low_ens_trend_mar_orig = trend_calc(siu_low_annavg_ens_mar,lat,lon,low_years); %2025-2(0,1)(9,0)0 MARCH
siu_high_ens_trend_mar_orig = trend_calc(siu_high_annavg_ens_mar,lat,lon,high_years); %19(5-6)0-2024 MARCH
siv_all_ens_trend_mar_orig = trend_calc(siv_all_annavg_ens_mar,lat,lon,all_years);
siv_low_ens_trend_sep_orig = trend_calc(siv_low_annavg_ens_sep,lat,lon,low_years);
siv_low_ens_trend_mar_orig = trend_calc(siv_low_annavg_ens_mar,lat,lon,low_years);
siv_high_ens_trend_mar_orig = trend_calc(siv_high_annavg_ens_mar,lat,lon,high_years);
u_rhs_geo_all_ens_trend_mar_orig = trend_calc(u_rhs_geo_all_annavg_ens_mar,lat,lon,all_years);
u_rhs_geo_low_ens_trend_sep_orig = trend_calc(u_rhs_geo_low_annavg_ens_sep,lat,lon,low_years);
u_rhs_geo_low_ens_trend_mar_orig = trend_calc(u_rhs_geo_low_annavg_ens_mar,lat,lon,low_years);
u_rhs_geo_high_ens_trend_mar_orig = trend_calc(u_rhs_geo_high_annavg_ens_mar,lat,lon,high_years);
u_rhs_ageo_all_ens_trend_mar_orig = trend_calc(u_rhs_ageo_all_annavg_ens_mar,lat,lon,all_years);
u_rhs_ageo_low_ens_trend_sep_orig = trend_calc(u_rhs_ageo_low_annavg_ens_sep,lat,lon,low_years);
u_rhs_ageo_low_ens_trend_mar_orig = trend_calc(u_rhs_ageo_low_annavg_ens_mar,lat,lon,low_years);
u_rhs_ageo_high_ens_trend_mar_orig = trend_calc(u_rhs_ageo_high_annavg_ens_mar,lat,lon,high_years);
u_rhs_total_all_ens_trend_mar_orig = trend_calc(u_rhs_total_all_annavg_ens_mar,lat,lon,all_years);
u_rhs_total_low_ens_trend_sep_orig = trend_calc(u_rhs_total_low_annavg_ens_sep,lat,lon,low_years);
u_rhs_total_low_ens_trend_mar_orig = trend_calc(u_rhs_total_low_annavg_ens_mar,lat,lon,low_years);
u_rhs_total_high_ens_trend_mar_orig = trend_calc(u_rhs_total_high_annavg_ens_mar,lat,lon,high_years);
v_rhs_geo_all_ens_trend_mar_orig = trend_calc(v_rhs_geo_all_annavg_ens_mar,lat,lon,all_years);
v_rhs_geo_low_ens_trend_sep_orig = trend_calc(v_rhs_geo_low_annavg_ens_sep,lat,lon,low_years);
v_rhs_geo_low_ens_trend_mar_orig = trend_calc(v_rhs_geo_low_annavg_ens_mar,lat,lon,low_years);
v_rhs_geo_high_ens_trend_mar_orig = trend_calc(v_rhs_geo_high_annavg_ens_mar,lat,lon,high_years);
v_rhs_ageo_all_ens_trend_mar_orig = trend_calc(v_rhs_ageo_all_annavg_ens_mar,lat,lon,all_years);
v_rhs_ageo_low_ens_trend_sep_orig = trend_calc(v_rhs_ageo_low_annavg_ens_sep,lat,lon,low_years);
v_rhs_ageo_low_ens_trend_mar_orig = trend_calc(v_rhs_ageo_low_annavg_ens_mar,lat,lon,low_years);
v_rhs_ageo_high_ens_trend_mar_orig = trend_calc(v_rhs_ageo_high_annavg_ens_mar,lat,lon,high_years);
v_rhs_total_all_ens_trend_mar_orig = trend_calc(v_rhs_total_all_annavg_ens_mar,lat,lon,all_years);
v_rhs_total_low_ens_trend_sep_orig = trend_calc(v_rhs_total_low_annavg_ens_sep,lat,lon,low_years);
v_rhs_total_low_ens_trend_mar_orig = trend_calc(v_rhs_total_low_annavg_ens_mar,lat,lon,low_years);
v_rhs_total_high_ens_trend_mar_orig = trend_calc(v_rhs_total_high_annavg_ens_mar,lat,lon,high_years);
sep_speed_low_annavg_ens_trend = trend_calc(sep_speed_low_annavg_ens,lat,lon,low_years);
mar_speed_low_annavg_ens_trend = trend_calc(mar_speed_low_annavg_ens,lat,lon,low_years);
u_ageo_from_siu_low_ens_trend_mar_orig = trend_calc(u_ageo_from_siu_low_annavg_ens_mar,lat,lon,low_years);
u_ageo_from_siu_low_ens_trend_sep_orig = trend_calc(u_ageo_from_siu_low_annavg_ens_sep,lat,lon,low_years);
v_ageo_from_siv_low_ens_trend_mar_orig = trend_calc(v_ageo_from_siv_low_annavg_ens_mar,lat,lon,low_years);
v_ageo_from_siv_low_ens_trend_sep_orig = trend_calc(v_ageo_from_siv_low_annavg_ens_sep,lat,lon,low_years);

% Angle correction: use angle (units of RADIANS) to find lat-lon grid
% directionality (east versus west; north versus south) from the original
% curvilinear grid. 
% Equations used: 
% x_variable_east = x_variable*cos(angle) - y_variable*sin(angle);
% y_variable_north = x_variable*sin(angle) + y_variable*cos(angle);
siu_all_ens_trend_mar = (siu_all_ens_trend_mar_orig.*cos(angle))-(siv_all_ens_trend_mar_orig.*sin(angle));
siu_high_ens_trend_mar = (siu_high_ens_trend_mar_orig.*cos(angle))-(siv_high_ens_trend_mar_orig.*sin(angle));
siu_low_ens_trend_sep = (siu_low_ens_trend_sep_orig.*cos(angle))-(siv_low_ens_trend_sep_orig.*sin(angle));
siu_low_ens_trend_mar = (siu_low_ens_trend_mar_orig.*cos(angle))-(siv_low_ens_trend_mar_orig.*sin(angle));
siv_all_ens_trend_mar = (siu_all_ens_trend_mar_orig.*sin(angle))+(siv_all_ens_trend_mar_orig.*cos(angle));
siv_high_ens_trend_mar = (siu_high_ens_trend_mar_orig.*sin(angle))+(siv_high_ens_trend_mar_orig.*cos(angle));
siv_low_ens_trend_sep = (siu_low_ens_trend_sep_orig.*sin(angle))+(siv_low_ens_trend_sep_orig.*cos(angle));
siv_low_ens_trend_mar = (siu_low_ens_trend_mar_orig.*sin(angle))+(siv_low_ens_trend_mar_orig.*cos(angle));
u_rhs_geo_all_ens_trend_mar = (u_rhs_geo_all_ens_trend_mar_orig.*cos(angle))-(v_rhs_geo_all_ens_trend_mar_orig.*sin(angle));
u_rhs_geo_high_ens_trend_mar = (u_rhs_geo_high_ens_trend_mar_orig.*cos(angle))-(v_rhs_geo_high_ens_trend_mar_orig.*sin(angle));
u_rhs_geo_low_ens_trend_sep = (u_rhs_geo_low_ens_trend_sep_orig.*cos(angle))-(v_rhs_geo_low_ens_trend_sep_orig.*sin(angle));
u_rhs_geo_low_ens_trend_mar = (u_rhs_geo_low_ens_trend_mar_orig.*cos(angle))-(v_rhs_geo_low_ens_trend_mar_orig.*sin(angle));
v_rhs_geo_all_ens_trend_mar = (u_rhs_geo_all_ens_trend_mar_orig.*sin(angle))+(v_rhs_geo_all_ens_trend_mar_orig.*cos(angle));
v_rhs_geo_high_ens_trend_mar = (u_rhs_geo_high_ens_trend_mar_orig.*sin(angle))+(v_rhs_geo_high_ens_trend_mar_orig.*cos(angle));
v_rhs_geo_low_ens_trend_sep = (u_rhs_geo_low_ens_trend_sep_orig.*sin(angle))+(v_rhs_geo_low_ens_trend_sep_orig.*cos(angle));
v_rhs_geo_low_ens_trend_mar = (u_rhs_geo_low_ens_trend_mar_orig.*sin(angle))+(v_rhs_geo_low_ens_trend_mar_orig.*cos(angle));
u_rhs_ageo_all_ens_trend_mar = (u_rhs_ageo_all_ens_trend_mar_orig.*cos(angle))-(v_rhs_ageo_all_ens_trend_mar_orig.*sin(angle));
u_rhs_ageo_high_ens_trend_mar = (u_rhs_ageo_high_ens_trend_mar_orig.*cos(angle))-(v_rhs_ageo_high_ens_trend_mar_orig.*sin(angle));
u_rhs_ageo_low_ens_trend_sep = (u_rhs_ageo_low_ens_trend_sep_orig.*cos(angle))-(v_rhs_ageo_low_ens_trend_sep_orig.*sin(angle));
u_rhs_ageo_low_ens_trend_mar = (u_rhs_ageo_low_ens_trend_mar_orig.*cos(angle))-(v_rhs_ageo_low_ens_trend_mar_orig.*sin(angle));
v_rhs_ageo_all_ens_trend_mar = (u_rhs_ageo_all_ens_trend_mar_orig.*sin(angle))+(v_rhs_ageo_all_ens_trend_mar_orig.*cos(angle));
v_rhs_ageo_high_ens_trend_mar = (u_rhs_ageo_high_ens_trend_mar_orig.*sin(angle))+(v_rhs_ageo_high_ens_trend_mar_orig.*cos(angle));
v_rhs_ageo_low_ens_trend_sep = (u_rhs_ageo_low_ens_trend_sep_orig.*sin(angle))+(v_rhs_ageo_low_ens_trend_sep_orig.*cos(angle));
v_rhs_ageo_low_ens_trend_mar = (u_rhs_ageo_low_ens_trend_mar_orig.*sin(angle))+(v_rhs_ageo_low_ens_trend_mar_orig.*cos(angle));
u_rhs_total_all_ens_trend_mar = (u_rhs_total_all_ens_trend_mar_orig.*cos(angle))-(v_rhs_total_all_ens_trend_mar_orig.*sin(angle));
u_rhs_total_high_ens_trend_mar = (u_rhs_total_high_ens_trend_mar_orig.*cos(angle))-(v_rhs_total_high_ens_trend_mar_orig.*sin(angle));
u_rhs_total_low_ens_trend_sep = (u_rhs_total_low_ens_trend_sep_orig.*cos(angle))-(v_rhs_total_low_ens_trend_sep_orig.*sin(angle));
u_rhs_total_low_ens_trend_mar = (u_rhs_total_low_ens_trend_mar_orig.*cos(angle))-(v_rhs_total_low_ens_trend_mar_orig.*sin(angle));
v_rhs_total_all_ens_trend_mar = (u_rhs_total_all_ens_trend_mar_orig.*sin(angle))+(v_rhs_total_all_ens_trend_mar_orig.*cos(angle));
v_rhs_total_high_ens_trend_mar = (u_rhs_total_high_ens_trend_mar_orig.*sin(angle))+(v_rhs_total_high_ens_trend_mar_orig.*cos(angle));
v_rhs_total_low_ens_trend_sep = (u_rhs_total_low_ens_trend_sep_orig.*sin(angle))+(v_rhs_total_low_ens_trend_sep_orig.*cos(angle));
v_rhs_total_low_ens_trend_mar = (u_rhs_total_low_ens_trend_mar_orig.*sin(angle))+(v_rhs_total_low_ens_trend_mar_orig.*cos(angle));
u_ageo_from_siu_low_ens_trend_mar = (u_ageo_from_siu_low_ens_trend_mar_orig.*cos(angle))-(v_ageo_from_siv_low_ens_trend_mar_orig.*sin(angle));
u_ageo_from_siu_low_ens_trend_sep = (u_ageo_from_siu_low_ens_trend_sep_orig.*cos(angle))-(v_ageo_from_siv_low_ens_trend_sep_orig.*sin(angle));
v_ageo_from_siv_low_ens_trend_mar = (u_ageo_from_siu_low_ens_trend_mar_orig.*sin(angle))+(v_ageo_from_siv_low_ens_trend_mar_orig.*cos(angle));
v_ageo_from_siv_low_ens_trend_sep = (u_ageo_from_siu_low_ens_trend_sep_orig.*sin(angle))+(v_ageo_from_siv_low_ens_trend_sep_orig.*cos(angle));

% Angle corrected trends converted from m s^-1 yr^-1 to km d^-1 dec^-1
conversion_factor_time = (1/1000) * (3600*24) * 10; %left to right: 1000 meters in 1km, 3600*24 seconds in one day, and 10 years in 1 decade.
siu_all_ens_trend_mar = siu_all_ens_trend_mar * conversion_factor_time;
siu_high_ens_trend_mar = siu_high_ens_trend_mar * conversion_factor_time;
siu_low_ens_trend_mar = siu_low_ens_trend_mar * conversion_factor_time;
siu_low_ens_trend_sep = siu_low_ens_trend_sep * conversion_factor_time;
siv_all_ens_trend_mar = siv_all_ens_trend_mar * conversion_factor_time;
siv_high_ens_trend_mar = siv_high_ens_trend_mar * conversion_factor_time;
siv_low_ens_trend_mar = siv_low_ens_trend_mar * conversion_factor_time;
siv_low_ens_trend_sep = siv_low_ens_trend_sep * conversion_factor_time;
u_rhs_geo_all_ens_trend_mar = u_rhs_geo_all_ens_trend_mar * conversion_factor_time;
u_rhs_geo_high_ens_trend_mar = u_rhs_geo_high_ens_trend_mar * conversion_factor_time;
u_rhs_geo_low_ens_trend_mar = u_rhs_geo_low_ens_trend_mar * conversion_factor_time;
u_rhs_geo_low_ens_trend_sep = u_rhs_geo_low_ens_trend_sep * conversion_factor_time;
v_rhs_geo_all_ens_trend_mar = v_rhs_geo_all_ens_trend_mar * conversion_factor_time;
v_rhs_geo_high_ens_trend_mar = v_rhs_geo_high_ens_trend_mar * conversion_factor_time;
v_rhs_geo_low_ens_trend_mar = v_rhs_geo_low_ens_trend_mar * conversion_factor_time;
v_rhs_geo_low_ens_trend_sep = v_rhs_geo_low_ens_trend_sep * conversion_factor_time;
u_rhs_ageo_all_ens_trend_mar = u_rhs_ageo_all_ens_trend_mar * conversion_factor_time;
u_rhs_ageo_high_ens_trend_mar = u_rhs_ageo_high_ens_trend_mar * conversion_factor_time;
u_rhs_ageo_low_ens_trend_mar = u_rhs_ageo_low_ens_trend_mar * conversion_factor_time;
u_rhs_ageo_low_ens_trend_sep = u_rhs_ageo_low_ens_trend_sep * conversion_factor_time;
v_rhs_ageo_all_ens_trend_mar = v_rhs_ageo_all_ens_trend_mar * conversion_factor_time;
v_rhs_ageo_high_ens_trend_mar = v_rhs_ageo_high_ens_trend_mar * conversion_factor_time;
v_rhs_ageo_low_ens_trend_mar = v_rhs_ageo_low_ens_trend_mar * conversion_factor_time;
v_rhs_ageo_low_ens_trend_sep = v_rhs_ageo_low_ens_trend_sep * conversion_factor_time;
u_rhs_total_all_ens_trend_mar = u_rhs_total_all_ens_trend_mar * conversion_factor_time;
u_rhs_total_high_ens_trend_mar = u_rhs_total_high_ens_trend_mar * conversion_factor_time;
u_rhs_total_low_ens_trend_mar = u_rhs_total_low_ens_trend_mar * conversion_factor_time;
u_rhs_total_low_ens_trend_sep = u_rhs_total_low_ens_trend_sep * conversion_factor_time;
v_rhs_total_all_ens_trend_mar = v_rhs_total_all_ens_trend_mar * conversion_factor_time;
v_rhs_total_high_ens_trend_mar = v_rhs_total_high_ens_trend_mar * conversion_factor_time;
v_rhs_total_low_ens_trend_mar = v_rhs_total_low_ens_trend_mar * conversion_factor_time;
v_rhs_total_low_ens_trend_sep = v_rhs_total_low_ens_trend_sep * conversion_factor_time;
sep_speed_low_annavg_ens_trend = sep_speed_low_annavg_ens_trend * conversion_factor_time;
mar_speed_low_annavg_ens_trend = mar_speed_low_annavg_ens_trend * conversion_factor_time;
u_ageo_from_siu_low_ens_trend_mar = u_ageo_from_siu_low_ens_trend_mar * conversion_factor_time;
u_ageo_from_siu_low_ens_trend_sep = u_ageo_from_siu_low_ens_trend_sep * conversion_factor_time;
v_ageo_from_siv_low_ens_trend_mar = v_ageo_from_siv_low_ens_trend_mar * conversion_factor_time;
v_ageo_from_siv_low_ens_trend_sep = v_ageo_from_siv_low_ens_trend_sep * conversion_factor_time;

% Angle-correct climatological data
siu_all_high_ens_clim_mar = (siu_all_high_ens_clim_orig_mar.*cos(angle))-(siv_all_high_ens_clim_orig_mar.*sin(angle));
siu_low_ens_clim_sep = (siu_low_ens_clim_orig_sep.*cos(angle))-(siv_low_ens_clim_orig_sep.*sin(angle));
siu_low_ens_clim_mar = (siu_low_ens_clim_orig_mar.*cos(angle))-(siv_low_ens_clim_orig_mar.*sin(angle));
siv_all_high_ens_clim_mar = (siu_all_high_ens_clim_orig_mar.*sin(angle))+(siv_all_high_ens_clim_orig_mar.*cos(angle));
siv_low_ens_clim_sep = (siu_low_ens_clim_orig_sep.*sin(angle))+(siv_low_ens_clim_orig_sep.*cos(angle));
siv_low_ens_clim_mar = (siu_low_ens_clim_orig_mar.*sin(angle))+(siv_low_ens_clim_orig_mar.*cos(angle));
u_rhs_geo_all_high_ens_clim_mar = (u_rhs_geo_all_high_ens_clim_orig_mar.*cos(angle))-(v_rhs_geo_all_high_ens_clim_orig_mar.*sin(angle));
u_rhs_geo_low_ens_clim_sep = (u_rhs_geo_low_ens_clim_orig_sep.*cos(angle))-(v_rhs_geo_low_ens_clim_orig_sep.*sin(angle));
u_rhs_geo_low_ens_clim_mar = (u_rhs_geo_low_ens_clim_orig_mar.*cos(angle))-(v_rhs_geo_low_ens_clim_orig_mar.*sin(angle));
u_rhs_ageo_all_high_ens_clim_mar = (u_rhs_ageo_all_high_ens_clim_orig_mar.*cos(angle))-(v_rhs_ageo_all_high_ens_clim_orig_mar.*sin(angle));
u_rhs_ageo_low_ens_clim_sep = (u_rhs_ageo_low_ens_clim_orig_sep.*cos(angle))-(v_rhs_ageo_low_ens_clim_orig_sep.*sin(angle));
u_rhs_ageo_low_ens_clim_mar = (u_rhs_ageo_low_ens_clim_orig_mar.*cos(angle))-(v_rhs_ageo_low_ens_clim_orig_mar.*sin(angle));
u_rhs_total_all_high_ens_clim_mar = (u_rhs_total_all_high_ens_clim_orig_mar.*cos(angle))-(v_rhs_total_all_high_ens_clim_orig_mar.*sin(angle));
u_rhs_total_low_ens_clim_sep = (u_rhs_total_low_ens_clim_orig_sep.*cos(angle))-(v_rhs_total_low_ens_clim_orig_sep.*sin(angle));
u_rhs_total_low_ens_clim_mar = (u_rhs_total_low_ens_clim_orig_mar.*cos(angle))-(v_rhs_total_low_ens_clim_orig_mar.*sin(angle));
v_rhs_geo_all_high_ens_clim_mar = (u_rhs_geo_all_high_ens_clim_orig_mar.*sin(angle))+(v_rhs_geo_all_high_ens_clim_orig_mar.*cos(angle));
v_rhs_geo_low_ens_clim_sep = (u_rhs_geo_low_ens_clim_orig_sep.*sin(angle))+(v_rhs_geo_low_ens_clim_orig_sep.*cos(angle));
v_rhs_geo_low_ens_clim_mar = (u_rhs_geo_low_ens_clim_orig_mar.*sin(angle))+(v_rhs_geo_low_ens_clim_orig_mar.*cos(angle));
v_rhs_ageo_all_high_ens_clim_mar = (u_rhs_ageo_all_high_ens_clim_orig_mar.*sin(angle))+(v_rhs_ageo_all_high_ens_clim_orig_mar.*cos(angle));
v_rhs_ageo_low_ens_clim_sep = (u_rhs_ageo_low_ens_clim_orig_sep.*sin(angle))+(v_rhs_ageo_low_ens_clim_orig_sep.*cos(angle));
v_rhs_ageo_low_ens_clim_mar = (u_rhs_ageo_low_ens_clim_orig_mar.*sin(angle))+(v_rhs_ageo_low_ens_clim_orig_mar.*cos(angle));
v_rhs_total_all_high_ens_clim_mar = (u_rhs_total_all_high_ens_clim_orig_mar.*sin(angle))+(v_rhs_total_all_high_ens_clim_orig_mar.*cos(angle));
v_rhs_total_low_ens_clim_sep = (u_rhs_total_low_ens_clim_orig_sep.*sin(angle))+(v_rhs_total_low_ens_clim_orig_sep.*cos(angle));
v_rhs_total_low_ens_clim_mar = (u_rhs_total_low_ens_clim_orig_mar.*sin(angle))+(v_rhs_total_low_ens_clim_orig_mar.*cos(angle));
u_ageo_from_siu_all_high_ens_clim_mar = (u_ageo_from_siu_all_high_ens_clim_orig_mar.*cos(angle))-(v_ageo_from_siv_all_high_ens_clim_orig_mar.*sin(angle));
u_ageo_from_siu_low_ens_clim_mar = (u_ageo_from_siu_low_ens_clim_orig_mar.*cos(angle))-(v_ageo_from_siv_low_ens_clim_orig_mar.*sin(angle));
u_ageo_from_siu_low_ens_clim_sep = (u_ageo_from_siu_low_ens_clim_orig_sep.*cos(angle))-(v_ageo_from_siv_low_ens_clim_orig_sep.*sin(angle));
v_ageo_from_siv_all_high_ens_clim_mar = (u_ageo_from_siu_all_high_ens_clim_orig_mar.*sin(angle))+(v_ageo_from_siv_all_high_ens_clim_orig_mar.*cos(angle));
v_ageo_from_siv_low_ens_clim_mar = (u_ageo_from_siu_low_ens_clim_orig_mar.*sin(angle))+(v_ageo_from_siv_low_ens_clim_orig_mar.*cos(angle));
v_ageo_from_siv_low_ens_clim_sep = (u_ageo_from_siu_low_ens_clim_orig_sep.*sin(angle))+(v_ageo_from_siv_low_ens_clim_orig_sep.*cos(angle));

% Climatologies: from m/s to km/day
conversion_factor_time_clim = conversion_factor_time/10; %No "per decade".
siu_all_high_ens_clim_mar = siu_all_high_ens_clim_mar * conversion_factor_time_clim;
siu_low_ens_clim_sep = siu_low_ens_clim_sep * conversion_factor_time_clim;
siu_low_ens_clim_mar = siu_low_ens_clim_mar * conversion_factor_time_clim;
siv_all_high_ens_clim_mar = siv_all_high_ens_clim_mar * conversion_factor_time_clim;
siv_low_ens_clim_sep = siv_low_ens_clim_sep * conversion_factor_time_clim;
siv_low_ens_clim_mar = siv_low_ens_clim_mar * conversion_factor_time_clim;
u_rhs_geo_all_high_ens_clim_mar = u_rhs_geo_all_high_ens_clim_mar * conversion_factor_time_clim;
u_rhs_geo_low_ens_clim_sep = u_rhs_geo_low_ens_clim_sep * conversion_factor_time_clim;
u_rhs_geo_low_ens_clim_mar = u_rhs_geo_low_ens_clim_mar * conversion_factor_time_clim;
v_rhs_geo_all_high_ens_clim_mar = v_rhs_geo_all_high_ens_clim_mar * conversion_factor_time_clim;
v_rhs_geo_low_ens_clim_sep = v_rhs_geo_low_ens_clim_sep * conversion_factor_time_clim;
v_rhs_geo_low_ens_clim_mar = v_rhs_geo_low_ens_clim_mar * conversion_factor_time_clim;
u_rhs_ageo_all_high_ens_clim_mar = u_rhs_ageo_all_high_ens_clim_mar * conversion_factor_time_clim;
u_rhs_ageo_low_ens_clim_sep = u_rhs_ageo_low_ens_clim_sep * conversion_factor_time_clim;
u_rhs_ageo_low_ens_clim_mar = u_rhs_ageo_low_ens_clim_mar * conversion_factor_time_clim;
v_rhs_ageo_all_high_ens_clim_mar = v_rhs_ageo_all_high_ens_clim_mar * conversion_factor_time_clim;
v_rhs_ageo_low_ens_clim_sep = v_rhs_ageo_low_ens_clim_sep * conversion_factor_time_clim;
v_rhs_ageo_low_ens_clim_mar = v_rhs_ageo_low_ens_clim_mar * conversion_factor_time_clim;
u_rhs_total_all_high_ens_clim_mar = u_rhs_total_all_high_ens_clim_mar * conversion_factor_time_clim;
u_rhs_total_low_ens_clim_sep = u_rhs_total_low_ens_clim_sep * conversion_factor_time_clim;
u_rhs_total_low_ens_clim_mar = u_rhs_total_low_ens_clim_mar * conversion_factor_time_clim;
v_rhs_total_all_high_ens_clim_mar = v_rhs_total_all_high_ens_clim_mar * conversion_factor_time_clim;
v_rhs_total_low_ens_clim_sep = v_rhs_total_low_ens_clim_sep * conversion_factor_time_clim;
v_rhs_total_low_ens_clim_mar = v_rhs_total_low_ens_clim_mar * conversion_factor_time_clim;
sep_speed_annavg_ens = sep_speed_annavg_ens * conversion_factor_time_clim;
mar_speed_annavg_ens = mar_speed_annavg_ens * conversion_factor_time_clim;
u_ageo_from_siu_all_high_ens_clim_mar = u_ageo_from_siu_all_high_ens_clim_mar * conversion_factor_time_clim;
u_ageo_from_siu_low_ens_clim_sep = u_ageo_from_siu_low_ens_clim_sep * conversion_factor_time_clim;
u_ageo_from_siu_low_ens_clim_mar = u_ageo_from_siu_low_ens_clim_mar * conversion_factor_time_clim;
v_ageo_from_siv_all_high_ens_clim_mar = v_ageo_from_siv_all_high_ens_clim_mar * conversion_factor_time_clim;
v_ageo_from_siv_low_ens_clim_sep = v_ageo_from_siv_low_ens_clim_sep * conversion_factor_time_clim;
v_ageo_from_siv_low_ens_clim_mar = v_ageo_from_siv_low_ens_clim_mar * conversion_factor_time_clim;

% Angle-corrected siu and siv data (NOT TRENDS):
siu_all_annavg_ens_sep_final = (siu_all_annavg_ens_sep(:,:,:).*cos(angle(:,:))) - (siv_all_annavg_ens_sep(:,:,:).*sin(angle(:,:)));
siu_all_annavg_ens_mar_final = (siu_all_annavg_ens_mar(:,:,:).*cos(angle(:,:))) - (siv_all_annavg_ens_mar(:,:,:).*sin(angle(:,:)));
siv_all_annavg_ens_sep_final = (siu_all_annavg_ens_sep(:,:,:).*sin(angle(:,:))) + (siv_all_annavg_ens_sep(:,:,:).*cos(angle(:,:)));
siv_all_annavg_ens_mar_final = (siu_all_annavg_ens_mar(:,:,:).*sin(angle(:,:))) + (siv_all_annavg_ens_mar(:,:,:).*cos(angle(:,:)));

% Convert siu and siv NON-trend data from m/s to km/d
siu_all_annavg_ens_sep_final = siu_all_annavg_ens_sep_final * conversion_factor_time_clim;
siu_all_annavg_ens_mar_final = siu_all_annavg_ens_mar_final * conversion_factor_time_clim;
siv_all_annavg_ens_sep_final = siv_all_annavg_ens_sep_final * conversion_factor_time_clim;
siv_all_annavg_ens_mar_final = siv_all_annavg_ens_mar_final * conversion_factor_time_clim;

% For plotting purposes, create lon_skipped, lat_skipped, and
% siu/siv_skipped variables (this way, vector fields won't saturate plot
% areas with arrows)
lon_skipped = lon(1:5:end,1:5:end);
lat_skipped = lat(1:5:end,1:5:end);
siu_all_ens_trend_mar_skipped = siu_all_ens_trend_mar(1:5:end,1:5:end);
siu_high_ens_trend_mar_skipped = siu_high_ens_trend_mar(1:5:end,1:5:end);
siu_low_ens_trend_sep_skipped = siu_low_ens_trend_sep(1:5:end,1:5:end);
siu_low_ens_trend_mar_skipped = siu_low_ens_trend_mar(1:5:end,1:5:end);
siv_all_ens_trend_skipped_mar = siv_all_ens_trend_mar(1:5:end,1:5:end);
siv_high_ens_trend_skipped_mar = siv_high_ens_trend_mar(1:5:end,1:5:end);
siv_low_ens_trend_skipped_sep = siv_low_ens_trend_sep(1:5:end,1:5:end);
siv_low_ens_trend_skipped_mar = siv_low_ens_trend_mar(1:5:end,1:5:end);

%% Ageostrophic component analysis: Atm+Ocn versus Stress
% Ageostrophic ensemble trends have been calculated from
% u_rhs_ageo_all_annavg_ens, etc. These have also been angle-corrected.
% Ageostrophic beginning time climatologies need to be calculated
% (1960-1980) and angle-corrected.
all_years = 1950:2100;
% chopped_flag = 1; %Like previous section, chopped flag==1 means year range will be 1960-2090. Else, year range is 1950-2100.

% Individual components (i.e., strairx(y)_d, strocnx(y)_d, and
% strintx(y)_d) have not been imported or angle corrected yet. 
in_strairx_var = [in_strairx_var_filename,'_div_fma']; % These variables account for steady state (i.e., divide by f*simass) on daily time scales.
in_strairy_var = [in_strairy_var_filename,'_div_fma']; % They were calculated in ensemble_force_divided_by_mass_times_coriolis.m
in_strocnx_var = [in_strocnx_var_filename,'_div_fma']; 
in_strocny_var = [in_strocny_var_filename,'_div_fma'];
in_strintx_var = [in_strintx_var_filename,'_div_fma'];
in_strinty_var = [in_strinty_var_filename,'_div_fma'];

% Define empty cell arrays for ageostrophic variables (by direction and
% season)
% September
u_strairy_rhs_all_sep_cell = cell(num_ens,1);
v_strairx_rhs_all_sep_cell = cell(num_ens,1);
u_strocny_rhs_all_sep_cell = cell(num_ens,1);
v_strocnx_rhs_all_sep_cell = cell(num_ens,1);
u_strinty_rhs_all_sep_cell = cell(num_ens,1);
v_strintx_rhs_all_sep_cell = cell(num_ens,1);
% March
u_strairy_rhs_all_mar_cell = cell(num_ens,1);
v_strairx_rhs_all_mar_cell = cell(num_ens,1);
u_strocny_rhs_all_mar_cell = cell(num_ens,1);
v_strocnx_rhs_all_mar_cell = cell(num_ens,1);
u_strinty_rhs_all_mar_cell = cell(num_ens,1);
v_strintx_rhs_all_mar_cell = cell(num_ens,1);

% Import atmosphere, ocean, and internal stress variables for each
% component direction.
% atm (x; strairy/fma in steady state)
for ens_index=1:num_ens
    % September
    u_strairy_rhs_squared_hist_annavg_sep = ncread([root_dir,in_strairy_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairy_var_filename,hist_file_end_sep],in_strairy_var);
    u_strairy_rhs_squared_fut_annavg_sep = ncread([root_dir,in_strairy_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairy_var_filename,fut_file_end_sep],in_strairy_var);
    u_strairy_rhs_all_sep_cell{ens_index,1}(:,:,:) = cat(3,u_strairy_rhs_squared_hist_annavg_sep,u_strairy_rhs_squared_fut_annavg_sep); %Concatenate historical and future data.
    clear u_strairy_rhs_squared_hist_annavg_sep
    clear u_strairy_rhs_squared_fut_annavg_sep
    % March
    u_strairy_rhs_squared_hist_annavg_mar = ncread([root_dir,in_strairy_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairy_var_filename,hist_file_end_mar],in_strairy_var);
    u_strairy_rhs_squared_fut_annavg_mar = ncread([root_dir,in_strairy_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairy_var_filename,fut_file_end_mar],in_strairy_var);
    u_strairy_rhs_all_mar_cell{ens_index,1}(:,:,:) = cat(3,u_strairy_rhs_squared_hist_annavg_mar,u_strairy_rhs_squared_fut_annavg_mar); %Concatenate historical and future data.
    clear u_strairy_rhs_squared_hist_annavg_mar
    clear u_strairy_rhs_squared_fut_annavg_mar    
end

% atm (y; strairx in steady state)
for ens_index=1:num_ens
    % September
    v_strairx_rhs_squared_hist_annavg_sep = ncread([root_dir,in_strairx_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairx_var_filename,hist_file_end_sep],in_strairx_var);
    v_strairx_rhs_squared_fut_annavg_sep = ncread([root_dir,in_strairx_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairx_var_filename,fut_file_end_sep],in_strairx_var);
    v_strairx_rhs_all_sep_cell{ens_index,1}(:,:,:) = cat(3,v_strairx_rhs_squared_hist_annavg_sep,v_strairx_rhs_squared_fut_annavg_sep); %Concatenate historical and future data.
    clear v_strairx_rhs_squared_hist_annavg_sep
    clear v_strairx_rhs_squared_fut_annavg_sep
    % March
    v_strairx_rhs_squared_hist_annavg_mar = ncread([root_dir,in_strairx_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairx_var_filename,hist_file_end_mar],in_strairx_var);
    v_strairx_rhs_squared_fut_annavg_mar = ncread([root_dir,in_strairx_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strairx_var_filename,fut_file_end_mar],in_strairx_var);
    v_strairx_rhs_all_mar_cell{ens_index,1}(:,:,:) = cat(3,v_strairx_rhs_squared_hist_annavg_mar,v_strairx_rhs_squared_fut_annavg_mar); %Concatenate historical and future data.
    clear v_strairx_rhs_squared_hist_annavg_mar
    clear v_strairx_rhs_squared_fut_annavg_mar
end

% ocn (x; strocny in steady state)
for ens_index=1:num_ens
    % September
    u_strocny_rhs_squared_hist_annavg_sep = ncread([root_dir,in_strocny_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocny_var_filename,hist_file_end_sep],in_strocny_var);
    u_strocny_rhs_squared_fut_annavg_sep = ncread([root_dir,in_strocny_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocny_var_filename,fut_file_end_sep],in_strocny_var);
    u_strocny_rhs_all_sep_cell{ens_index,1}(:,:,:) = cat(3,u_strocny_rhs_squared_hist_annavg_sep,u_strocny_rhs_squared_fut_annavg_sep); %Concatenate historical and future data.
    clear u_strocny_rhs_squared_hist_annavg_sep
    clear u_strocny_rhs_squared_fut_annavg_sep
    % March
    u_strocny_rhs_squared_hist_annavg_mar = ncread([root_dir,in_strocny_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocny_var_filename,hist_file_end_mar],in_strocny_var);
    u_strocny_rhs_squared_fut_annavg_mar = ncread([root_dir,in_strocny_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocny_var_filename,fut_file_end_mar],in_strocny_var);
    u_strocny_rhs_all_mar_cell{ens_index,1}(:,:,:) = cat(3,u_strocny_rhs_squared_hist_annavg_mar,u_strocny_rhs_squared_fut_annavg_mar); %Concatenate historical and future data.
    clear u_strocny_rhs_squared_hist_annavg_mar
    clear u_strocny_rhs_squared_fut_annavg_mar
end

% ocn (y; strocnx in steady state)
for ens_index=1:num_ens
    % September 
    v_strocnx_rhs_squared_hist_annavg_sep = ncread([root_dir,in_strocnx_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocnx_var_filename,hist_file_end_sep],in_strocnx_var);
    v_strocnx_rhs_squared_fut_annavg_sep = ncread([root_dir,in_strocnx_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocnx_var_filename,fut_file_end_sep],in_strocnx_var);
    v_strocnx_rhs_all_sep_cell{ens_index,1}(:,:,:) = cat(3,v_strocnx_rhs_squared_hist_annavg_sep,v_strocnx_rhs_squared_fut_annavg_sep); %Concatenate historical and future data.
    clear v_strocnx_rhs_squared_hist_annavg_sep
    clear v_strocnx_rhs_squared_fut_annavg_sep
    % March
    v_strocnx_rhs_squared_hist_annavg_mar = ncread([root_dir,in_strocnx_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocnx_var_filename,hist_file_end_mar],in_strocnx_var);
    v_strocnx_rhs_squared_fut_annavg_mar = ncread([root_dir,in_strocnx_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strocnx_var_filename,fut_file_end_mar],in_strocnx_var);
    v_strocnx_rhs_all_mar_cell{ens_index,1}(:,:,:) = cat(3,v_strocnx_rhs_squared_hist_annavg_mar,v_strocnx_rhs_squared_fut_annavg_mar); %Concatenate historical and future data.
    clear v_strocnx_rhs_squared_hist_annavg_mar
    clear v_strocnx_rhs_squared_fut_annavg_mar
end

% stress (x; strinty in steady state)
for ens_index=1:num_ens
    % September
    u_strinty_rhs_squared_hist_annavg_sep = ncread([root_dir,in_strinty_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strinty_var_filename,hist_file_end_sep],in_strinty_var);
    u_strinty_rhs_squared_fut_annavg_sep = ncread([root_dir,in_strinty_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strinty_var_filename,fut_file_end_sep],in_strinty_var);
    u_strinty_rhs_all_sep_cell{ens_index,1}(:,:,:) = cat(3,u_strinty_rhs_squared_hist_annavg_sep,u_strinty_rhs_squared_fut_annavg_sep); % Concatenate historical and future data.
    clear u_strinty_rhs_squared_hist_annavg_sep
    clear u_strinty_rhs_squared_fut_annavg_sep
    % March
    u_strinty_rhs_squared_hist_annavg_mar = ncread([root_dir,in_strinty_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strinty_var_filename,hist_file_end_mar],in_strinty_var);
    u_strinty_rhs_squared_fut_annavg_mar = ncread([root_dir,in_strinty_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strinty_var_filename,fut_file_end_mar],in_strinty_var);
    u_strinty_rhs_all_mar_cell{ens_index,1}(:,:,:) = cat(3,u_strinty_rhs_squared_hist_annavg_mar,u_strinty_rhs_squared_fut_annavg_mar); % Concatenate historical and future data.
    clear u_strinty_rhs_squared_hist_annavg_mar
    clear u_strinty_rhs_squared_fut_annavg_mar
end

% stress (y; strintx in steady state)
for ens_index=1:num_ens
    % September
    v_strintx_rhs_squared_hist_annavg_sep = ncread([root_dir,in_strintx_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strintx_var_filename,hist_file_end_sep],in_strintx_var);
    v_strintx_rhs_squared_fut_annavg_sep = ncread([root_dir,in_strintx_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strintx_var_filename,fut_file_end_sep],in_strintx_var);
    v_strintx_rhs_all_sep_cell{ens_index,1}(:,:,:) = cat(3,v_strintx_rhs_squared_hist_annavg_sep,v_strintx_rhs_squared_fut_annavg_sep); % Concatenate historical and future data.
    clear v_strintx_rhs_squared_hist_annavg_sep
    clear v_strintx_rhs_squared_fut_annavg_sep
    % March
    v_strintx_rhs_squared_hist_annavg_mar = ncread([root_dir,in_strintx_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strintx_var_filename,hist_file_end_mar],in_strintx_var);
    v_strintx_rhs_squared_fut_annavg_mar = ncread([root_dir,in_strintx_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_strintx_var_filename,fut_file_end_mar],in_strintx_var);
    v_strintx_rhs_all_mar_cell{ens_index,1}(:,:,:) = cat(3,v_strintx_rhs_squared_hist_annavg_mar,v_strintx_rhs_squared_fut_annavg_mar); % Concatenate historical and future data.
    clear v_strintx_rhs_squared_hist_annavg_mar
    clear v_strintx_rhs_squared_fut_annavg_mar
end

% Now, filter out cells based on dist2coast (150km or NaN) and lat/lon designations
for index=1:num_ens %Finally, cycle through each ensemble member.
    for row_idx=1:num_lons %Then, cycle through each longitude.
        for col_idx=1:num_lats %First, cycle through each latitude.
            if dist2coast(row_idx,col_idx)<150 || isnan(dist2coast(row_idx,col_idx)) %If the grid in question is within 150km of land or dist2coast at that grid is a NaN value.
                u_strairy_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strairy_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
            end
        end
    end
end
% Final filter: NaN out grids outside of the "Tandon domain"
% (doi:10.1029/2017JC013697)
for index=1:num_ens %Finally, cycle through each ensemble member.
    for row_idx=1:num_lons %Then, cycle through each longitude.
        for col_idx=1:num_lats %First, cycle through each latitude.
            if (lon(row_idx,col_idx)>=103 && lon(row_idx,col_idx)<=236) && lat(row_idx,col_idx)<68 %Latitudes should be >=68. If not, nan out the corresponding row_idx,col_idx value in speed.
                u_strairy_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN; 
                u_strairy_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
            elseif (lon(row_idx,col_idx)>236 || lon(row_idx,col_idx)<103) && lat(row_idx,col_idx)<79 %Latitudes in this longitude range should be at least 79 degN
                u_strairy_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN; 
                u_strairy_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
            elseif isnan(lon(row_idx,col_idx)) || isnan(lat(row_idx,col_idx)) % If lat or lon is a NaN value, the corresponding ageostrophic variables must also be NaN'ed.
                u_strairy_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strairy_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strairx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strocny_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                v_strocnx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
                u_strinty_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;           
                v_strintx_rhs_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
            end
        end
    end
end

% Steady state equations: du/dt = -(g/f)(dH/dy) + u_strairy + u_strocny +
% u_strinty (ageo_u = u_strairy + u_strocny + u_strinty). dv/dt =
% (g/f)(dH/dx) - v_strairx - v_strocnx - v_strintx (ageo_v = -v_strairx -
% v_strocnx - v_strintx). Here, we calculate atm + ocn for each direction
% component. In the case of v, we multiply all ageostrophic components by
% negative 1 within their cell arrays.
% Define empty cell arrays to hold atm + ocn sums for each ensemble member.
u_strairy_strocny_rhs_all_sep_cell = cell(num_ens,1);
v_strairx_strocnx_rhs_all_sep_cell = cell(num_ens,1);
u_ind_ageo_rhs_all_sep_cell = cell(num_ens,1);
v_ind_ageo_rhs_all_sep_cell = cell(num_ens,1);
u_strairy_strocny_rhs_all_mar_cell = cell(num_ens,1);
v_strairx_strocnx_rhs_all_mar_cell = cell(num_ens,1);
u_ind_ageo_rhs_all_mar_cell = cell(num_ens,1);
v_ind_ageo_rhs_all_mar_cell = cell(num_ens,1);

for ens_index=1:num_ens %Calculate atm+ocn and multiply all v components by negative 1 for each ensemble member.
    u_strairy_strocny_rhs_all_sep_cell{ens_index,1}(:,:,:) = u_strairy_rhs_all_sep_cell{ens_index,1}(:,:,:) + u_strocny_rhs_all_sep_cell{ens_index,1}(:,:,:);
    % Sum of u steady state ageostrophic components. NOT USED IN FINAL
    % CALCULATIONS FOR PAPER.
    u_ind_ageo_rhs_all_sep_cell{ens_index,1}(:,:,:) = u_strairy_rhs_all_sep_cell{ens_index,1}(:,:,:) + u_strocny_rhs_all_sep_cell{ens_index,1}(:,:,:) + u_strinty_rhs_all_sep_cell{ens_index,1}(:,:,:);
    v_strairx_rhs_all_sep_cell{ens_index,1}(:,:,:) = -1 * v_strairx_rhs_all_sep_cell{ens_index,1}(:,:,:);
    v_strocnx_rhs_all_sep_cell{ens_index,1}(:,:,:) = -1 * v_strocnx_rhs_all_sep_cell{ens_index,1}(:,:,:);
    v_strintx_rhs_all_sep_cell{ens_index,1}(:,:,:) = -1 * v_strintx_rhs_all_sep_cell{ens_index,1}(:,:,:);
    v_strairx_strocnx_rhs_all_sep_cell{ens_index,1}(:,:,:) = v_strairx_rhs_all_sep_cell{ens_index,1}(:,:,:) + v_strocnx_rhs_all_sep_cell{ens_index,1}(:,:,:);
    % Sum of v steady state ageostrophic components. NOT USED IN FINAL
    % CALCULATIONS FOR PAPER.
    v_ind_ageo_rhs_all_sep_cell{ens_index,1}(:,:,:) = v_strairx_rhs_all_sep_cell{ens_index,1}(:,:,:) + v_strocnx_rhs_all_sep_cell{ens_index,1}(:,:,:) + v_strintx_rhs_all_sep_cell{ens_index,1}(:,:,:);
    u_strairy_strocny_rhs_all_mar_cell{ens_index,1}(:,:,:) = u_strairy_rhs_all_mar_cell{ens_index,1}(:,:,:) + u_strocny_rhs_all_mar_cell{ens_index,1}(:,:,:);
    % Sum of v steady state ageostrophic components. NOT USED IN FINAL
    % CALCULATIONS FOR PAPER.
    u_ind_ageo_rhs_all_mar_cell{ens_index,1}(:,:,:) = u_strairy_rhs_all_mar_cell{ens_index,1}(:,:,:) + u_strocny_rhs_all_mar_cell{ens_index,1}(:,:,:) + u_strinty_rhs_all_mar_cell{ens_index,1}(:,:,:);
    v_strairx_rhs_all_mar_cell{ens_index,1}(:,:,:) = -1 * v_strairx_rhs_all_mar_cell{ens_index,1}(:,:,:);
    v_strocnx_rhs_all_mar_cell{ens_index,1}(:,:,:) = -1 * v_strocnx_rhs_all_mar_cell{ens_index,1}(:,:,:);
    v_strintx_rhs_all_mar_cell{ens_index,1}(:,:,:) = -1 * v_strintx_rhs_all_mar_cell{ens_index,1}(:,:,:);
    v_strairx_strocnx_rhs_all_mar_cell{ens_index,1}(:,:,:) = v_strairx_rhs_all_mar_cell{ens_index,1}(:,:,:) + v_strocnx_rhs_all_mar_cell{ens_index,1}(:,:,:);
    % Sum of v steady state ageostrophic components. NOT USED IN FINAL
    % CALCULATIONS FOR PAPER.
    v_ind_ageo_rhs_all_mar_cell{ens_index,1}(:,:,:) = v_strairx_rhs_all_mar_cell{ens_index,1}(:,:,:) + v_strocnx_rhs_all_mar_cell{ens_index,1}(:,:,:) + v_strintx_rhs_all_mar_cell{ens_index,1}(:,:,:);
end

% Ensemble average calculation time! 
% Define empty 3D NaN arrays, where the length of all_years is the third
% dimension.
u_strairy_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
v_strairx_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
u_strocny_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
v_strocnx_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
u_strinty_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
v_strintx_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
u_strairy_strocny_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
v_strairx_strocnx_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
u_ind_ageo_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
v_ind_ageo_all_annavg_sep_ens = nan(num_lons,num_lats,length(all_years));
u_strairy_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
v_strairx_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
u_strocny_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
v_strocnx_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
u_strinty_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
v_strintx_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
u_strairy_strocny_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
v_strairx_strocnx_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
u_ind_ageo_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));
v_ind_ageo_all_annavg_mar_ens = nan(num_lons,num_lats,length(all_years));

% Calculate ensemble average for each year.
for index=1:length(all_years)
    % Define temporary arrays to hold spatial data for each ensemble (for a
    % given year="index").
    temp_u_strairy_sep = nan(num_lons,num_lats,num_ens);
    temp_v_strairx_sep = nan(num_lons,num_lats,num_ens);
    temp_u_strocny_sep = nan(num_lons,num_lats,num_ens);
    temp_v_strocnx_sep = nan(num_lons,num_lats,num_ens);
    temp_u_strinty_sep = nan(num_lons,num_lats,num_ens);
    temp_v_strintx_sep = nan(num_lons,num_lats,num_ens);
    temp_u_strairy_strocny_sep = nan(num_lons,num_lats,num_ens);
    temp_v_strairx_strocnx_sep = nan(num_lons,num_lats,num_ens);
    temp_u_ind_ageo_sep = nan(num_lons,num_lats,num_ens);
    temp_v_ind_ageo_sep = nan(num_lons,num_lats,num_ens);
    temp_u_strairy_mar = nan(num_lons,num_lats,num_ens);
    temp_v_strairx_mar = nan(num_lons,num_lats,num_ens);
    temp_u_strocny_mar = nan(num_lons,num_lats,num_ens);
    temp_v_strocnx_mar = nan(num_lons,num_lats,num_ens);
    temp_u_strinty_mar = nan(num_lons,num_lats,num_ens);
    temp_v_strintx_mar = nan(num_lons,num_lats,num_ens);
    temp_u_strairy_strocny_mar = nan(num_lons,num_lats,num_ens);
    temp_v_strairx_strocnx_mar = nan(num_lons,num_lats,num_ens);
    temp_u_ind_ageo_mar = nan(num_lons,num_lats,num_ens);
    temp_v_ind_ageo_mar = nan(num_lons,num_lats,num_ens);
    % Extract "index" time slice from each cell, place it in
    % temp_*(:,:,cell_index)
    for cell_index=1:num_ens
        temp_u_strairy_sep(:,:,cell_index) = u_strairy_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_v_strairx_sep(:,:,cell_index) = v_strairx_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_u_strocny_sep(:,:,cell_index) = u_strocny_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_v_strocnx_sep(:,:,cell_index) = v_strocnx_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_u_strinty_sep(:,:,cell_index) = u_strinty_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_v_strintx_sep(:,:,cell_index) = v_strintx_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_u_strairy_strocny_sep(:,:,cell_index) = u_strairy_strocny_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_v_strairx_strocnx_sep(:,:,cell_index) = v_strairx_strocnx_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_u_ind_ageo_sep(:,:,cell_index) = u_ind_ageo_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_v_ind_ageo_sep(:,:,cell_index) = v_ind_ageo_rhs_all_sep_cell{cell_index,1}(:,:,index);
        temp_u_strairy_mar(:,:,cell_index) = u_strairy_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_v_strairx_mar(:,:,cell_index) = v_strairx_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_u_strocny_mar(:,:,cell_index) = u_strocny_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_v_strocnx_mar(:,:,cell_index) = v_strocnx_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_u_strinty_mar(:,:,cell_index) = u_strinty_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_v_strintx_mar(:,:,cell_index) = v_strintx_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_u_strairy_strocny_mar(:,:,cell_index) = u_strairy_strocny_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_v_strairx_strocnx_mar(:,:,cell_index) = v_strairx_strocnx_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_u_ind_ageo_mar(:,:,cell_index) = u_ind_ageo_rhs_all_mar_cell{cell_index,1}(:,:,index);
        temp_v_ind_ageo_mar(:,:,cell_index) = v_ind_ageo_rhs_all_mar_cell{cell_index,1}(:,:,index);
    end
    % Calculate the ensemble average by taking the mean of temp_* along the
    % time dimension. Place the resulting average in *_ens(:,:,index).
    % NaN's are not included in the average calculation.
    u_strairy_all_annavg_sep_ens(:,:,index) = mean(temp_u_strairy_sep,3,'omitnan');
    v_strairx_all_annavg_sep_ens(:,:,index) = mean(temp_v_strairx_sep,3,'omitnan');
    u_strocny_all_annavg_sep_ens(:,:,index) = mean(temp_u_strocny_sep,3,'omitnan');
    v_strocnx_all_annavg_sep_ens(:,:,index) = mean(temp_v_strocnx_sep,3,'omitnan');
    u_strinty_all_annavg_sep_ens(:,:,index) = mean(temp_u_strinty_sep,3,'omitnan');
    v_strintx_all_annavg_sep_ens(:,:,index) = mean(temp_v_strintx_sep,3,'omitnan');
    u_strairy_strocny_all_annavg_sep_ens(:,:,index) = mean(temp_u_strairy_strocny_sep,3,'omitnan');
    v_strairx_strocnx_all_annavg_sep_ens(:,:,index) = mean(temp_v_strairx_strocnx_sep,3,'omitnan');
    u_ind_ageo_all_annavg_sep_ens(:,:,index) = mean(temp_u_ind_ageo_sep,3,'omitnan');
    v_ind_ageo_all_annavg_sep_ens(:,:,index) = mean(temp_v_ind_ageo_sep,3,'omitnan');
    u_strairy_all_annavg_mar_ens(:,:,index) = mean(temp_u_strairy_mar,3,'omitnan');
    v_strairx_all_annavg_mar_ens(:,:,index) = mean(temp_v_strairx_mar,3,'omitnan');
    u_strocny_all_annavg_mar_ens(:,:,index) = mean(temp_u_strocny_mar,3,'omitnan');
    v_strocnx_all_annavg_mar_ens(:,:,index) = mean(temp_v_strocnx_mar,3,'omitnan');
    u_strinty_all_annavg_mar_ens(:,:,index) = mean(temp_u_strinty_mar,3,'omitnan');
    v_strintx_all_annavg_mar_ens(:,:,index) = mean(temp_v_strintx_mar,3,'omitnan');
    u_strairy_strocny_all_annavg_mar_ens(:,:,index) = mean(temp_u_strairy_strocny_mar,3,'omitnan');
    v_strairx_strocnx_all_annavg_mar_ens(:,:,index) = mean(temp_v_strairx_strocnx_mar,3,'omitnan');
    u_ind_ageo_all_annavg_mar_ens(:,:,index) = mean(temp_u_ind_ageo_mar,3,'omitnan');
    v_ind_ageo_all_annavg_mar_ens(:,:,index) = mean(temp_v_ind_ageo_mar,3,'omitnan');
end

% Climatologies: 1950-1970 (1-21 for All, high years) and 2025-2035 (66-86 for Low years) for u and v direction ageostrophic individual terms.
% Total ageostrophic climatologies were calculated and angle-corrected 2 sections above: variable names are
% u_rhs_ageo_all_high_ens_clim, v_rhs_ageo_all_high_ens_clim,
% u_rhs_ageo_low_ens_clim, and v_rhs_ageo_low_ens_clim.
% By-term climatologies are NOT angle-corrected.
u_strairy_all_high_ens_clim_orig_mar = mean(u_strairy_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
u_strairy_low_ens_clim_orig_sep = mean(u_strairy_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
u_strairy_low_ens_clim_orig_mar = mean(u_strairy_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
v_strairx_all_high_ens_clim_orig_mar = mean(v_strairx_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
v_strairx_low_ens_clim_orig_sep = mean(v_strairx_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
v_strairx_low_ens_clim_orig_mar = mean(v_strairx_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
u_strocny_all_high_ens_clim_orig_mar = mean(u_strocny_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
u_strocny_low_ens_clim_orig_sep = mean(u_strocny_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
u_strocny_low_ens_clim_orig_mar = mean(u_strocny_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
v_strocnx_all_high_ens_clim_orig_mar = mean(v_strocnx_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
v_strocnx_low_ens_clim_orig_sep = mean(v_strocnx_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
v_strocnx_low_ens_clim_orig_mar = mean(v_strocnx_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
u_strinty_all_high_ens_clim_orig_mar = mean(u_strinty_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
u_strinty_low_ens_clim_orig_sep = mean(u_strinty_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
u_strinty_low_ens_clim_orig_mar = mean(u_strinty_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
v_strintx_all_high_ens_clim_orig_mar = mean(v_strintx_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
v_strintx_low_ens_clim_orig_sep = mean(v_strintx_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
v_strintx_low_ens_clim_orig_mar = mean(v_strintx_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
% atm+ocn sum climatologies
u_strairy_strocny_all_high_ens_clim_orig_mar = mean(u_strairy_strocny_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
u_strairy_strocny_low_ens_clim_orig_sep = mean(u_strairy_strocny_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
u_strairy_strocny_low_ens_clim_orig_mar = mean(u_strairy_strocny_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');
v_strairx_strocnx_all_high_ens_clim_orig_mar = mean(v_strairx_strocnx_all_annavg_mar_ens(:,:,high_years_indices),3,'omitnan');
v_strairx_strocnx_low_ens_clim_orig_sep = mean(v_strairx_strocnx_all_annavg_sep_ens(:,:,low_years_indices),3,'omitnan');
v_strairx_strocnx_low_ens_clim_orig_mar = mean(v_strairx_strocnx_all_annavg_mar_ens(:,:,low_years_indices),3,'omitnan');

if chopped_flag==1 %Only looking at 1960-2090, not 1950-2100.
    % Now, for each cell element, cut off the first and last 10 years. 
    % Cut off first and last 10 years from each component model
    % % Cut off first 10 and last 10 years
    all_years = 1960:2090; % Redefine all_years from 1:131 in previous section so subsetting in this section works.
    all_years = all_years-1949; %First value in all_years is now 11 (neglecting 1-10 to cut them off)

    % Chop off first and last 10 years.
    u_strairy_all_annavg_sep_ens = u_strairy_all_annavg_sep_ens(:,:,all_years);
    v_strairx_all_annavg_sep_ens = v_strairx_all_annavg_sep_ens(:,:,all_years);
    u_strocny_all_annavg_sep_ens = u_strocny_all_annavg_sep_ens(:,:,all_years);
    v_strocnx_all_annavg_sep_ens = v_strocnx_all_annavg_sep_ens(:,:,all_years);
    u_strinty_all_annavg_sep_ens = u_strinty_all_annavg_sep_ens(:,:,all_years);
    v_strintx_all_annavg_sep_ens = v_strintx_all_annavg_sep_ens(:,:,all_years);
    u_strairy_strocny_all_annavg_sep_ens = u_strairy_strocny_all_annavg_sep_ens(:,:,all_years);
    v_strairx_strocnx_all_annavg_sep_ens = v_strairx_strocnx_all_annavg_sep_ens(:,:,all_years);
    u_ind_ageo_all_annavg_sep_ens = u_ind_ageo_all_annavg_sep_ens(:,:,all_years);
    v_ind_ageo_all_annavg_sep_ens = v_ind_ageo_all_annavg_sep_ens(:,:,all_years);
    u_strairy_all_annavg_mar_ens = u_strairy_all_annavg_mar_ens(:,:,all_years);
    v_strairx_all_annavg_mar_ens = v_strairx_all_annavg_mar_ens(:,:,all_years);
    u_strocny_all_annavg_mar_ens = u_strocny_all_annavg_mar_ens(:,:,all_years);
    v_strocnx_all_annavg_mar_ens = v_strocnx_all_annavg_mar_ens(:,:,all_years);
    u_strinty_all_annavg_mar_ens = u_strinty_all_annavg_mar_ens(:,:,all_years);
    v_strintx_all_annavg_mar_ens = v_strintx_all_annavg_mar_ens(:,:,all_years);
    u_strairy_strocny_all_annavg_mar_ens = u_strairy_strocny_all_annavg_mar_ens(:,:,all_years);
    v_strairx_strocnx_all_annavg_mar_ens = v_strairx_strocnx_all_annavg_mar_ens(:,:,all_years);
    u_ind_ageo_all_annavg_mar_ens = u_ind_ageo_all_annavg_mar_ens(:,:,all_years);
    v_ind_ageo_all_annavg_mar_ens = v_ind_ageo_all_annavg_mar_ens(:,:,all_years);

    % Change all_years so all_years(1)=1 and corresponds to 1960.
    all_years = all_years-10;

    % Define high (increasing speed) and low (decreasing drift speed)
    % years.
    high_years = 1960:2024; 
    high_years = high_years-1959; %The first index==1 now.
    low_years = 2025:2090; 
    low_years = low_years-1959;
else %Not getting rid of any ageostrophic data, so just redefining all_years and specifying when drift speed increases (high_years) or decreases (low_years)
    all_years = all_years-1949;
    high_years = 1950:2024;
    high_years = high_years-1949;
    low_years = 2025:2100;
    low_years = low_years-1949;
end

% Extract appropriate data for each year chunk and each ensemble variable:
% atmospheric stress in the u direction for the steady state assumption.
u_strairy_low_annavg_ens_sep = u_strairy_all_annavg_sep_ens(:,:,low_years); %SEP, low_years (projected)
u_strairy_low_annavg_ens_mar = u_strairy_all_annavg_mar_ens(:,:,low_years); %MARCH, low_years (projected)
u_strairy_high_annavg_ens_mar = u_strairy_all_annavg_mar_ens(:,:,high_years); %MARCH, high_years (historical)
% atmospheric stress in the v direction.
v_strairx_low_annavg_ens_sep = v_strairx_all_annavg_sep_ens(:,:,low_years);
v_strairx_low_annavg_ens_mar = v_strairx_all_annavg_mar_ens(:,:,low_years);
v_strairx_high_annavg_ens_mar = v_strairx_all_annavg_mar_ens(:,:,high_years);
% oceanic stress in the u direction
u_strocny_low_annavg_ens_sep = u_strocny_all_annavg_sep_ens(:,:,low_years);
u_strocny_low_annavg_ens_mar = u_strocny_all_annavg_mar_ens(:,:,low_years);
u_strocny_high_annavg_ens_mar = u_strocny_all_annavg_mar_ens(:,:,high_years);
% oceanic stress in the v direction
v_strocnx_low_annavg_ens_sep = v_strocnx_all_annavg_sep_ens(:,:,low_years);
v_strocnx_low_annavg_ens_mar = v_strocnx_all_annavg_mar_ens(:,:,low_years);
v_strocnx_high_annavg_ens_mar = v_strocnx_all_annavg_mar_ens(:,:,high_years);
% internal stress in the u direction
u_strinty_low_annavg_ens_sep = u_strinty_all_annavg_sep_ens(:,:,low_years);
u_strinty_low_annavg_ens_mar = u_strinty_all_annavg_mar_ens(:,:,low_years);
u_strinty_high_annavg_ens_mar = u_strinty_all_annavg_mar_ens(:,:,high_years);
% internal stress in the v direction
v_strintx_low_annavg_ens_sep = v_strintx_all_annavg_sep_ens(:,:,low_years);
v_strintx_low_annavg_ens_mar = v_strintx_all_annavg_mar_ens(:,:,low_years);
v_strintx_high_annavg_ens_mar = v_strintx_all_annavg_mar_ens(:,:,high_years);
% atm+ocn stresses in the u direction
u_strairy_strocny_low_annavg_ens_sep = u_strairy_strocny_all_annavg_sep_ens(:,:,low_years);
u_strairy_strocny_low_annavg_ens_mar = u_strairy_strocny_all_annavg_mar_ens(:,:,low_years);
u_strairy_strocny_high_annavg_ens_mar = u_strairy_strocny_all_annavg_mar_ens(:,:,high_years);
% atm+ocn stresses in the v direction
v_strairx_strocnx_low_annavg_ens_sep = v_strairx_strocnx_all_annavg_sep_ens(:,:,low_years);
v_strairx_strocnx_low_annavg_ens_mar = v_strairx_strocnx_all_annavg_mar_ens(:,:,low_years);
v_strairx_strocnx_high_annavg_ens_mar = v_strairx_strocnx_all_annavg_mar_ens(:,:,high_years);

% ageostrophic sums of individual component variables (NOT USED IN FINAL
% CALCULATIONS FOR THE PAPER).
u_ind_ageo_low_annavg_ens_sep = u_ind_ageo_all_annavg_sep_ens(:,:,low_years);
u_ind_ageo_low_annavg_ens_mar = u_ind_ageo_all_annavg_mar_ens(:,:,low_years);
u_ind_ageo_high_annavg_ens_mar = u_ind_ageo_all_annavg_mar_ens(:,:,high_years);
v_ind_ageo_low_annavg_ens_sep = v_ind_ageo_all_annavg_sep_ens(:,:,low_years);
v_ind_ageo_low_annavg_ens_mar = v_ind_ageo_all_annavg_mar_ens(:,:,low_years);
v_ind_ageo_high_annavg_ens_mar = v_ind_ageo_all_annavg_mar_ens(:,:,high_years);

% In September, some lat/lon combinations may only have a couple of non-ice
% free years. We delete/NaN out these grid points to avoid erroneous trend
% calculations; pixels with only 5 years of non-NaN
% entries will be completely NaN'ed before trends are calculated.

for row_idx=1:num_lons %Finally, cycle by longitude.
    for col_idx=1:num_lats %Then, cycle by latitude.
        num_non_nans = 0; %Initialize a counter for non-NaN entries.
        for time_index=1:length(low_years) %First, cycle through each time step in "low_years"
            if ~isnan(siu_low_annavg_ens_sep(row_idx,col_idx,time_index)) %If the siu_low_annavg_... value at this lat/lon for a year=time_index is non-NaN, add to the counter num_non_nans
                num_non_nans = num_non_nans + 1;
            end
        end
        if num_non_nans<=5 && col_idx>5 && row_idx>5%If there are only 5 or less non-NaN years, we NaN these out now.
            for time_index=1:length(low_years) %Iterate through each year and NaN out non-NaN values
                if ~isnan(siu_low_annavg_ens_sep(row_idx,col_idx,time_index))
                    u_strairy_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_strairx_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_strocny_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_strocnx_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_strinty_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_strintx_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    u_strairy_strocny_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                    v_strairx_strocnx_low_annavg_ens_sep(row_idx,col_idx,time_index) = NaN;
                end
            end
        end
    end
end


% Calculate trends:
% Transpose low and high years for ease of calculation: 
high_years = high_years';
low_years = low_years';
all_years = all_years';

% Calculate trends using function trend_calc
u_strairy_all_ens_trend_orig_mar = trend_calc(u_strairy_all_annavg_mar_ens,lat,lon,all_years);
u_strairy_low_ens_trend_orig_sep = trend_calc(u_strairy_low_annavg_ens_sep,lat,lon,low_years);
u_strairy_low_ens_trend_orig_mar = trend_calc(u_strairy_low_annavg_ens_mar,lat,lon,low_years);
u_strairy_high_ens_trend_orig_mar = trend_calc(u_strairy_high_annavg_ens_mar,lat,lon,high_years);
v_strairx_all_ens_trend_orig_mar = trend_calc(v_strairx_all_annavg_mar_ens,lat,lon,all_years);
v_strairx_low_ens_trend_orig_sep = trend_calc(v_strairx_low_annavg_ens_sep,lat,lon,low_years);
v_strairx_low_ens_trend_orig_mar = trend_calc(v_strairx_low_annavg_ens_mar,lat,lon,low_years);
v_strairx_high_ens_trend_orig_mar = trend_calc(v_strairx_high_annavg_ens_mar,lat,lon,high_years);
u_strocny_all_ens_trend_orig_mar = trend_calc(u_strocny_all_annavg_mar_ens,lat,lon,all_years);
u_strocny_low_ens_trend_orig_sep = trend_calc(u_strocny_low_annavg_ens_sep,lat,lon,low_years);
u_strocny_low_ens_trend_orig_mar = trend_calc(u_strocny_low_annavg_ens_mar,lat,lon,low_years);
u_strocny_high_ens_trend_orig_mar = trend_calc(u_strocny_high_annavg_ens_mar,lat,lon,high_years);
v_strocnx_all_ens_trend_orig_mar = trend_calc(v_strocnx_all_annavg_mar_ens,lat,lon,all_years);
v_strocnx_low_ens_trend_orig_sep = trend_calc(v_strocnx_low_annavg_ens_sep,lat,lon,low_years);
v_strocnx_low_ens_trend_orig_mar = trend_calc(v_strocnx_low_annavg_ens_mar,lat,lon,low_years);
v_strocnx_high_ens_trend_orig_mar = trend_calc(v_strocnx_high_annavg_ens_mar,lat,lon,high_years);
u_strinty_all_ens_trend_orig_mar = trend_calc(u_strinty_all_annavg_mar_ens,lat,lon,all_years);
u_strinty_low_ens_trend_orig_sep = trend_calc(u_strinty_low_annavg_ens_sep,lat,lon,low_years);
u_strinty_low_ens_trend_orig_mar = trend_calc(u_strinty_low_annavg_ens_mar,lat,lon,low_years);
u_strinty_high_ens_trend_orig_mar = trend_calc(u_strinty_high_annavg_ens_mar,lat,lon,high_years);
v_strintx_all_ens_trend_orig_mar = trend_calc(v_strintx_all_annavg_mar_ens,lat,lon,all_years); 
v_strintx_low_ens_trend_orig_sep = trend_calc(v_strintx_low_annavg_ens_sep,lat,lon,low_years);
v_strintx_low_ens_trend_orig_mar = trend_calc(v_strintx_low_annavg_ens_mar,lat,lon,low_years);
v_strintx_high_ens_trend_orig_mar = trend_calc(v_strintx_high_annavg_ens_mar,lat,lon,high_years);
u_strairy_strocny_all_ens_trend_orig_mar = trend_calc(u_strairy_strocny_all_annavg_mar_ens,lat,lon,all_years);
u_strairy_strocny_low_ens_trend_orig_sep = trend_calc(u_strairy_strocny_low_annavg_ens_sep,lat,lon,low_years);
u_strairy_strocny_low_ens_trend_orig_mar = trend_calc(u_strairy_strocny_low_annavg_ens_mar,lat,lon,low_years);
u_strairy_strocny_high_ens_trend_orig_mar = trend_calc(u_strairy_strocny_high_annavg_ens_mar,lat,lon,high_years);
v_strairx_strocnx_all_ens_trend_orig_mar = trend_calc(v_strairx_strocnx_all_annavg_mar_ens,lat,lon,all_years);
v_strairx_strocnx_low_ens_trend_orig_sep = trend_calc(v_strairx_strocnx_low_annavg_ens_sep,lat,lon,low_years);
v_strairx_strocnx_low_ens_trend_orig_mar = trend_calc(v_strairx_strocnx_low_annavg_ens_mar,lat,lon,low_years);
v_strairx_strocnx_high_ens_trend_orig_mar = trend_calc(v_strairx_strocnx_high_annavg_ens_mar,lat,lon,high_years);
u_ind_ageo_all_ens_trend_orig_mar = trend_calc(u_ind_ageo_all_annavg_mar_ens,lat,lon,all_years);
u_ind_ageo_low_ens_trend_orig_sep = trend_calc(u_ind_ageo_low_annavg_ens_sep,lat,lon,low_years);
u_ind_ageo_low_ens_trend_orig_mar = trend_calc(u_ind_ageo_low_annavg_ens_mar,lat,lon,low_years);
u_ind_ageo_high_ens_trend_orig_mar = trend_calc(u_ind_ageo_high_annavg_ens_mar,lat,lon,high_years);
v_ind_ageo_all_ens_trend_orig_mar = trend_calc(v_ind_ageo_all_annavg_mar_ens,lat,lon,all_years);
v_ind_ageo_low_ens_trend_orig_sep = trend_calc(v_ind_ageo_low_annavg_ens_sep,lat,lon,low_years);
v_ind_ageo_low_ens_trend_orig_mar = trend_calc(v_ind_ageo_low_annavg_ens_mar,lat,lon,low_years);
v_ind_ageo_high_ens_trend_orig_mar = trend_calc(v_ind_ageo_high_annavg_ens_mar,lat,lon,high_years);

% Angle correct: new_u = old_u*cos(angle)-old_v*sin(angle) and new_v =
% old_u*sin(angle)+old_v*cos(angle)
u_strairy_all_ens_trend_mar = (u_strairy_all_ens_trend_orig_mar.*cos(angle))-(v_strairx_all_ens_trend_orig_mar.*sin(angle));
u_strairy_low_ens_trend_sep = (u_strairy_low_ens_trend_orig_sep.*cos(angle))-(v_strairx_low_ens_trend_orig_sep.*sin(angle));
u_strairy_low_ens_trend_mar = (u_strairy_low_ens_trend_orig_mar.*cos(angle))-(v_strairx_low_ens_trend_orig_mar.*sin(angle));
u_strairy_high_ens_trend_mar = (u_strairy_high_ens_trend_orig_mar.*cos(angle))-(v_strairx_high_ens_trend_orig_mar.*sin(angle));
v_strairx_all_ens_trend_mar = (u_strairy_all_ens_trend_orig_mar.*sin(angle))+(v_strairx_all_ens_trend_orig_mar.*cos(angle));
v_strairx_low_ens_trend_sep = (u_strairy_low_ens_trend_orig_sep.*sin(angle))+(v_strairx_low_ens_trend_orig_sep.*cos(angle));
v_strairx_low_ens_trend_mar = (u_strairy_low_ens_trend_orig_mar.*sin(angle))+(v_strairx_low_ens_trend_orig_mar.*cos(angle));
v_strairx_high_ens_trend_mar = (u_strairy_high_ens_trend_orig_mar.*sin(angle))+(v_strairx_high_ens_trend_orig_mar.*cos(angle));
u_strocny_all_ens_trend_mar = (u_strocny_all_ens_trend_orig_mar.*cos(angle))-(v_strocnx_all_ens_trend_orig_mar.*sin(angle));
u_strocny_low_ens_trend_sep = (u_strocny_low_ens_trend_orig_sep.*cos(angle))-(v_strocnx_low_ens_trend_orig_sep.*sin(angle));
u_strocny_low_ens_trend_mar = (u_strocny_low_ens_trend_orig_mar.*cos(angle))-(v_strocnx_low_ens_trend_orig_mar.*sin(angle));
u_strocny_high_ens_trend_mar = (u_strocny_high_ens_trend_orig_mar.*cos(angle))-(v_strocnx_high_ens_trend_orig_mar.*sin(angle));
v_strocnx_all_ens_trend_mar = (u_strocny_all_ens_trend_orig_mar.*sin(angle))+(v_strocnx_all_ens_trend_orig_mar.*cos(angle));
v_strocnx_low_ens_trend_sep = (u_strocny_low_ens_trend_orig_sep.*sin(angle))+(v_strocnx_low_ens_trend_orig_sep.*cos(angle));
v_strocnx_low_ens_trend_mar = (u_strocny_low_ens_trend_orig_mar.*sin(angle))+(v_strocnx_low_ens_trend_orig_mar.*cos(angle));
v_strocnx_high_ens_trend_mar = (u_strocny_high_ens_trend_orig_mar.*sin(angle))+(v_strocnx_high_ens_trend_orig_mar.*cos(angle));
u_strinty_all_ens_trend_mar = (u_strinty_all_ens_trend_orig_mar.*cos(angle))-(v_strintx_all_ens_trend_orig_mar.*sin(angle));
u_strinty_low_ens_trend_sep = (u_strinty_low_ens_trend_orig_sep.*cos(angle))-(v_strintx_low_ens_trend_orig_sep.*sin(angle));
u_strinty_low_ens_trend_mar = (u_strinty_low_ens_trend_orig_mar.*cos(angle))-(v_strintx_low_ens_trend_orig_mar.*sin(angle));
u_strinty_high_ens_trend_mar = (u_strinty_high_ens_trend_orig_mar.*cos(angle))-(v_strintx_high_ens_trend_orig_mar.*sin(angle));
v_strintx_all_ens_trend_mar = (u_strinty_all_ens_trend_orig_mar.*sin(angle))+(v_strintx_all_ens_trend_orig_mar.*cos(angle));
v_strintx_low_ens_trend_sep = (u_strinty_low_ens_trend_orig_sep.*sin(angle))+(v_strintx_low_ens_trend_orig_sep.*cos(angle));
v_strintx_low_ens_trend_mar = (u_strinty_low_ens_trend_orig_mar.*sin(angle))+(v_strintx_low_ens_trend_orig_mar.*cos(angle));
v_strintx_high_ens_trend_mar = (u_strinty_high_ens_trend_orig_mar.*sin(angle))+(v_strintx_high_ens_trend_orig_mar.*cos(angle));
u_strairy_plus_strocny_all_ens_trend_mar = (u_strairy_strocny_all_ens_trend_orig_mar.*cos(angle))-(v_strairx_strocnx_all_ens_trend_orig_mar.*sin(angle));
u_strairy_plus_strocny_low_ens_trend_sep = (u_strairy_strocny_low_ens_trend_orig_sep.*cos(angle))-(v_strairx_strocnx_low_ens_trend_orig_sep.*sin(angle));
u_strairy_plus_strocny_low_ens_trend_mar = (u_strairy_strocny_low_ens_trend_orig_mar.*cos(angle))-(v_strairx_strocnx_low_ens_trend_orig_mar.*sin(angle));
u_strairy_plus_strocny_high_ens_trend_mar = (u_strairy_strocny_high_ens_trend_orig_mar.*cos(angle))-(v_strairx_strocnx_high_ens_trend_orig_mar.*sin(angle));
v_strairx_plus_strocnx_all_ens_trend_mar = (u_strairy_strocny_all_ens_trend_orig_mar.*sin(angle))+(v_strairx_strocnx_all_ens_trend_orig_mar.*cos(angle));
v_strairx_plus_strocnx_low_ens_trend_sep = (u_strairy_strocny_low_ens_trend_orig_sep.*sin(angle))+(v_strairx_strocnx_low_ens_trend_orig_sep.*cos(angle));
v_strairx_plus_strocnx_low_ens_trend_mar = (u_strairy_strocny_low_ens_trend_orig_mar.*sin(angle))+(v_strairx_strocnx_low_ens_trend_orig_mar.*cos(angle));
v_strairx_plus_strocnx_high_ens_trend_mar = (u_strairy_strocny_high_ens_trend_orig_mar.*sin(angle))+(v_strairx_strocnx_high_ens_trend_orig_mar.*cos(angle));
u_ind_ageo_all_ens_trend_mar = (u_ind_ageo_all_ens_trend_orig_mar.*cos(angle)) - (v_ind_ageo_all_ens_trend_orig_mar.*sin(angle));
u_ind_ageo_low_ens_trend_sep = (u_ind_ageo_low_ens_trend_orig_sep.*cos(angle)) - (v_ind_ageo_low_ens_trend_orig_sep.*sin(angle));
u_ind_ageo_low_ens_trend_mar = (u_ind_ageo_low_ens_trend_orig_mar.*cos(angle)) - (v_ind_ageo_low_ens_trend_orig_mar.*sin(angle));
u_ind_ageo_high_ens_trend_mar = (u_ind_ageo_high_ens_trend_orig_mar.*cos(angle)) - (v_ind_ageo_high_ens_trend_orig_mar.*sin(angle));
v_ind_ageo_all_ens_trend_mar = (u_ind_ageo_all_ens_trend_orig_mar.*sin(angle)) + (v_ind_ageo_all_ens_trend_orig_mar.*cos(angle));
v_ind_ageo_low_ens_trend_sep = (u_ind_ageo_low_ens_trend_orig_sep.*sin(angle)) + (v_ind_ageo_low_ens_trend_orig_sep.*cos(angle));
v_ind_ageo_low_ens_trend_mar = (u_ind_ageo_low_ens_trend_orig_mar.*sin(angle)) + (v_ind_ageo_low_ens_trend_orig_mar.*cos(angle));
v_ind_ageo_high_ens_trend_mar = (u_ind_ageo_high_ens_trend_orig_mar.*sin(angle)) + (v_ind_ageo_high_ens_trend_orig_mar.*cos(angle));
% u_ageo_ens_clim = (u_ageo_ens_clim_all_high_orig.*cos(angle))-(v_ageo_ens_clim_all_high_orig.*sin(angle));
% v_ageo_ens_clim = (u_ageo_ens_clim_all_high_orig.*sin(angle))+(v_ageo_ens_clim_all_high_orig.*cos(angle));

% Conversion: from m/(s*yr) to km d^-1 dec^-1. conversion_factor_time is
% defined in the previous section.
u_strairy_all_ens_trend_mar = u_strairy_all_ens_trend_mar * conversion_factor_time;
u_strairy_low_ens_trend_sep = u_strairy_low_ens_trend_sep * conversion_factor_time;
u_strairy_low_ens_trend_mar = u_strairy_low_ens_trend_mar * conversion_factor_time;
u_strairy_high_ens_trend_mar = u_strairy_high_ens_trend_mar * conversion_factor_time;
v_strairx_all_ens_trend_mar = v_strairx_all_ens_trend_mar * conversion_factor_time;
v_strairx_low_ens_trend_sep = v_strairx_low_ens_trend_sep * conversion_factor_time;
v_strairx_low_ens_trend_mar = v_strairx_low_ens_trend_mar * conversion_factor_time;
v_strairx_high_ens_trend_mar = v_strairx_high_ens_trend_mar * conversion_factor_time;
u_strocny_all_ens_trend_mar = u_strocny_all_ens_trend_mar * conversion_factor_time;
u_strocny_low_ens_trend_sep = u_strocny_low_ens_trend_sep * conversion_factor_time;
u_strocny_low_ens_trend_mar = u_strocny_low_ens_trend_mar * conversion_factor_time;
u_strocny_high_ens_trend_mar = u_strocny_high_ens_trend_mar * conversion_factor_time;
v_strocnx_all_ens_trend_mar = v_strocnx_all_ens_trend_mar * conversion_factor_time;
v_strocnx_low_ens_trend_sep = v_strocnx_low_ens_trend_sep * conversion_factor_time;
v_strocnx_low_ens_trend_mar = v_strocnx_low_ens_trend_mar * conversion_factor_time;
v_strocnx_high_ens_trend_mar = v_strocnx_high_ens_trend_mar * conversion_factor_time;
u_strinty_all_ens_trend_mar = u_strinty_all_ens_trend_mar * conversion_factor_time;
u_strinty_low_ens_trend_sep = u_strinty_low_ens_trend_sep * conversion_factor_time;
u_strinty_low_ens_trend_mar = u_strinty_low_ens_trend_mar * conversion_factor_time;
u_strinty_high_ens_trend_mar = u_strinty_high_ens_trend_mar * conversion_factor_time;
v_strintx_all_ens_trend_mar = v_strintx_all_ens_trend_mar * conversion_factor_time;
v_strintx_low_ens_trend_sep = v_strintx_low_ens_trend_sep * conversion_factor_time;
v_strintx_low_ens_trend_mar = v_strintx_low_ens_trend_mar * conversion_factor_time;
v_strintx_high_ens_trend_mar = v_strintx_high_ens_trend_mar * conversion_factor_time;
u_strairy_plus_strocny_all_ens_trend_mar = u_strairy_plus_strocny_all_ens_trend_mar * conversion_factor_time;
u_strairy_plus_strocny_low_ens_trend_sep = u_strairy_plus_strocny_low_ens_trend_sep * conversion_factor_time;
u_strairy_plus_strocny_low_ens_trend_mar = u_strairy_plus_strocny_low_ens_trend_mar * conversion_factor_time;
u_strairy_plus_strocny_high_ens_trend_mar = u_strairy_plus_strocny_high_ens_trend_mar * conversion_factor_time;
v_strairx_plus_strocnx_all_ens_trend_mar = v_strairx_plus_strocnx_all_ens_trend_mar * conversion_factor_time;
v_strairx_plus_strocnx_low_ens_trend_sep = v_strairx_plus_strocnx_low_ens_trend_sep * conversion_factor_time;
v_strairx_plus_strocnx_low_ens_trend_mar = v_strairx_plus_strocnx_low_ens_trend_mar * conversion_factor_time;
v_strairx_plus_strocnx_high_ens_trend_mar = v_strairx_plus_strocnx_high_ens_trend_mar * conversion_factor_time;

%Climatologies: 
% Angle correct
u_strairy_all_high_ens_clim_mar = (u_strairy_all_high_ens_clim_orig_mar.*cos(angle)) - (v_strairx_all_high_ens_clim_orig_mar.*sin(angle));
u_strairy_low_ens_clim_sep = (u_strairy_low_ens_clim_orig_sep.*cos(angle)) - (v_strairx_low_ens_clim_orig_sep.*sin(angle));
u_strairy_low_ens_clim_mar = (u_strairy_low_ens_clim_orig_mar.*cos(angle)) - (v_strairx_low_ens_clim_orig_mar.*sin(angle));
v_strairx_all_high_ens_clim_mar = (u_strairy_all_high_ens_clim_orig_mar.*sin(angle)) + (v_strairx_all_high_ens_clim_orig_mar.*cos(angle));
v_strairx_low_ens_clim_sep = (u_strairy_low_ens_clim_orig_sep.*sin(angle)) + (v_strairx_low_ens_clim_orig_sep.*cos(angle));
v_strairx_low_ens_clim_mar = (u_strairy_low_ens_clim_orig_mar.*sin(angle)) + (v_strairx_low_ens_clim_orig_mar.*cos(angle));
u_strocny_all_high_ens_clim_mar = (u_strocny_all_high_ens_clim_orig_mar.*cos(angle)) - (v_strocnx_all_high_ens_clim_orig_mar.*sin(angle));
u_strocny_low_ens_clim_sep = (u_strocny_low_ens_clim_orig_sep.*cos(angle)) - (v_strocnx_low_ens_clim_orig_sep.*sin(angle));
u_strocny_low_ens_clim_mar = (u_strocny_low_ens_clim_orig_mar.*cos(angle)) - (v_strocnx_low_ens_clim_orig_mar.*sin(angle));
v_strocnx_all_high_ens_clim_mar = (u_strocny_all_high_ens_clim_orig_mar.*sin(angle)) + (v_strocnx_all_high_ens_clim_orig_mar.*cos(angle));
v_strocnx_low_ens_clim_sep = (u_strocny_low_ens_clim_orig_sep.*sin(angle)) + (v_strocnx_low_ens_clim_orig_sep.*cos(angle));
v_strocnx_low_ens_clim_mar = (u_strocny_low_ens_clim_orig_mar.*sin(angle)) + (v_strocnx_low_ens_clim_orig_mar.*cos(angle));
u_strinty_all_high_ens_clim_mar = (u_strinty_all_high_ens_clim_orig_mar.*cos(angle)) - (v_strintx_all_high_ens_clim_orig_mar.*sin(angle));
u_strinty_low_ens_clim_sep = (u_strinty_low_ens_clim_orig_sep.*cos(angle)) - (v_strintx_low_ens_clim_orig_sep.*sin(angle));
u_strinty_low_ens_clim_mar = (u_strinty_low_ens_clim_orig_mar.*cos(angle)) - (v_strintx_low_ens_clim_orig_mar.*sin(angle));
v_strintx_all_high_ens_clim_mar = (u_strinty_all_high_ens_clim_orig_mar.*sin(angle)) + (v_strintx_all_high_ens_clim_orig_mar.*cos(angle));
v_strintx_low_ens_clim_sep = (u_strinty_low_ens_clim_orig_sep.*sin(angle)) + (v_strintx_low_ens_clim_orig_sep.*cos(angle));
v_strintx_low_ens_clim_mar = (u_strinty_low_ens_clim_orig_mar.*sin(angle)) + (v_strintx_low_ens_clim_orig_mar.*cos(angle));
u_strairy_strocny_all_high_ens_clim_mar = (u_strairy_strocny_all_high_ens_clim_orig_mar.*cos(angle)) - (v_strairx_strocnx_all_high_ens_clim_orig_mar.*sin(angle));
u_strairy_strocny_low_ens_clim_sep = (u_strairy_strocny_low_ens_clim_orig_sep.*cos(angle)) - (v_strairx_strocnx_low_ens_clim_orig_sep.*sin(angle));
u_strairy_strocny_low_ens_clim_mar = (u_strairy_strocny_low_ens_clim_orig_mar.*cos(angle)) - (v_strairx_strocnx_low_ens_clim_orig_mar.*sin(angle));
v_strairx_strocnx_all_high_ens_clim_mar = (u_strairy_strocny_all_high_ens_clim_orig_mar.*sin(angle)) + (v_strairx_strocnx_all_high_ens_clim_orig_mar.*cos(angle));
v_strairx_strocnx_low_ens_clim_sep = (u_strairy_strocny_low_ens_clim_orig_sep.*sin(angle)) + (v_strairx_strocnx_low_ens_clim_orig_sep.*cos(angle));
v_strairx_strocnx_low_ens_clim_mar = (u_strairy_strocny_low_ens_clim_orig_mar.*sin(angle)) + (v_strairx_strocnx_low_ens_clim_orig_mar.*cos(angle));

% Conversion for climatologies (m/s to km/d). conversion_factor_time_clim
% is defined in the previous section.
u_strairy_all_high_ens_clim_mar = u_strairy_all_high_ens_clim_mar * conversion_factor_time_clim;
u_strairy_low_ens_clim_sep = u_strairy_low_ens_clim_sep * conversion_factor_time_clim;
u_strairy_low_ens_clim_mar = u_strairy_low_ens_clim_mar * conversion_factor_time_clim;
v_strairx_all_high_ens_clim_mar = v_strairx_all_high_ens_clim_mar * conversion_factor_time_clim;
v_strairx_low_ens_clim_sep = v_strairx_low_ens_clim_sep * conversion_factor_time_clim;
v_strairx_low_ens_clim_mar = v_strairx_low_ens_clim_mar * conversion_factor_time_clim;
u_strocny_all_high_ens_clim_mar = u_strocny_all_high_ens_clim_mar * conversion_factor_time_clim;
u_strocny_low_ens_clim_sep = u_strocny_low_ens_clim_sep * conversion_factor_time_clim;
u_strocny_low_ens_clim_mar = u_strocny_low_ens_clim_mar * conversion_factor_time_clim;
v_strocnx_all_high_ens_clim_mar = v_strocnx_all_high_ens_clim_mar * conversion_factor_time_clim;
v_strocnx_low_ens_clim_sep = v_strocnx_low_ens_clim_sep * conversion_factor_time_clim;
v_strocnx_low_ens_clim_mar = v_strocnx_low_ens_clim_mar * conversion_factor_time_clim;
u_strinty_all_high_ens_clim_mar = u_strinty_all_high_ens_clim_mar * conversion_factor_time_clim;
u_strinty_low_ens_clim_sep = u_strinty_low_ens_clim_sep * conversion_factor_time_clim;
u_strinty_low_ens_clim_mar = u_strinty_low_ens_clim_mar * conversion_factor_time_clim;
v_strintx_all_high_ens_clim_mar = v_strintx_all_high_ens_clim_mar * conversion_factor_time_clim;
v_strintx_low_ens_clim_sep = v_strintx_low_ens_clim_sep * conversion_factor_time_clim;
v_strintx_low_ens_clim_mar = v_strintx_low_ens_clim_mar * conversion_factor_time_clim;
u_strairy_strocny_all_high_ens_clim_mar = u_strairy_strocny_all_high_ens_clim_mar * conversion_factor_time_clim;
u_strairy_strocny_low_ens_clim_sep = u_strairy_strocny_low_ens_clim_sep * conversion_factor_time_clim;
u_strairy_strocny_low_ens_clim_mar = u_strairy_strocny_low_ens_clim_mar * conversion_factor_time_clim;
v_strairx_strocnx_all_high_ens_clim_mar = v_strairx_strocnx_all_high_ens_clim_mar * conversion_factor_time_clim;
v_strairx_strocnx_low_ens_clim_sep = v_strairx_strocnx_low_ens_clim_sep * conversion_factor_time_clim;
v_strairx_strocnx_low_ens_clim_mar = v_strairx_strocnx_low_ens_clim_mar * conversion_factor_time_clim;

%% SSH versus tilt (i.e., geostrophic force)
% Tilt force is associated with changes in sea surface height that cause
% sea ice to move down slope (from higher sea surface heights to lower sea
% surface heights). It's equivalent to the atmosphere's pressure gradient
% force. Here, we're comparing sea surface height trends to tilt force
% trends in each direction. Tilt force has already been angle corrected and
% SSH (a scalar field) does not need to undergo this process. 
% Reset all_years
all_years = 1950:2100;
chopped_flag = 0;

in_tilt_force_mag_var = 'tilt_force_mag';

% Create empty cell arrays
ssh_all_sep_cell = cell(num_ens,1);
ssh_all_mar_cell = cell(num_ens,1);
% tilt_mag_all_cell = cell(num_ens,1);

% SSH
for ens_index=1:num_ens
    ssh_hist_annavg_sep = ncread([root_dir,in_ssh_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),pop_file_middle,in_ssh_var,hist_file_end_sep],in_ssh_var);
    ssh_fut_annavg_sep = ncread([root_dir,in_ssh_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),pop_file_middle,in_ssh_var,fut_file_end_sep],in_ssh_var);
    ssh_all_sep_cell{ens_index,1}(:,:,:) = cat(3,ssh_hist_annavg_sep,ssh_fut_annavg_sep);
    %ssh_all_sep_cell{ens_index,1}(:,:,:) =
    %ssh_all_sep_cell{ens_index,1}(:,:,:) * (1/100); % SSH in the original
    %files is in units of cm. If we want meters, this is how we convert.
    clear ssh_hist_annavg_sep
    clear ssh_fut_annavg_sep
    ssh_hist_annavg_mar = ncread([root_dir,in_ssh_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),pop_file_middle,in_ssh_var,hist_file_end_mar],in_ssh_var);
    ssh_fut_annavg_mar = ncread([root_dir,in_ssh_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),pop_file_middle,in_ssh_var,fut_file_end_mar],in_ssh_var);
    ssh_all_mar_cell{ens_index,1}(:,:,:) = cat(3,ssh_hist_annavg_mar,ssh_fut_annavg_mar);
    %ssh_all_mar_cell{ens_index,1}(:,:,:) =
    %ssh_all_mar_cell{ens_index,1}(:,:,:) * (1/100); % SSH in the original
    %files is in units of cm. If we want meters, this is how we convert.
    clear ssh_hist_annavg_mar
    clear ssh_fut_annavg_mar
end

% % Import Tilt magnitude (calculated from daily tilt, equal to the square root of
% % each component squared, summed).
% for ens_index=1:num_ens
%     tilt_mag_hist_annavg = ncread([root_dir,in_ssh_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),pop_file_middle,in_ssh_var,hist_file_end],in_tilt_force_mag_var);
%     tilt_mag_fut_annavg = ncread([root_dir,in_ssh_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),pop_file_middle,in_ssh_var,fut_file_end],in_tilt_force_mag_var);
%     tilt_mag_all_cell{ens_index,1}(:,:,:) = cat(3,tilt_mag_hist_annavg,tilt_mag_fut_annavg);
% end

% Now, filter out cells based on dist2coast and lat/lon designations
for index=1:num_ens
    for row_idx=1:num_lons
        for col_idx=1:num_lats
            if dist2coast(row_idx,col_idx)<150 || isnan(dist2coast(row_idx,col_idx))
                ssh_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                ssh_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
%                 tilt_mag_all_cell{index,1}(row_idx,col_idx,:)=NaN;
            end
        end
    end
end

% Final filter: NaN out grids outside of the "Tandon domain"
for index=1:num_ens
    for row_idx=1:num_lons
        for col_idx=1:num_lats
            if (lon(row_idx,col_idx)>=103 && lon(row_idx,col_idx)<=236) && lat(row_idx,col_idx)<68 %Latitudes should be >=68. If not, nan out the corresponding row_idx,col_idx value in speed.
                ssh_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                ssh_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
%                 tilt_mag_all_cell{index,1}(row_idx,col_idx,:)=NaN;
            elseif (lon(row_idx,col_idx)>236 || lon(row_idx,col_idx)<103) && lat(row_idx,col_idx)<79 %Latitudes in this longitude range should be at least 79 degN
                ssh_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                ssh_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
%                 tilt_mag_all_cell{index,1}(row_idx,col_idx,:)=NaN;
            elseif isnan(lon(row_idx,col_idx)) || isnan(lat(row_idx,col_idx))
                ssh_all_sep_cell{index,1}(row_idx,col_idx,:)=NaN;
                ssh_all_mar_cell{index,1}(row_idx,col_idx,:)=NaN;
%                 tilt_mag_all_cell{index,1}(row_idx,col_idx,:)=NaN;
            end
        end
    end
end

% Ensemble average calculation time! 
ssh_all_annavg_ens_sep = nan(num_lons,num_lats,length(all_years));
ssh_all_annavg_ens_mar = nan(num_lons,num_lats,length(all_years));
% tilt_mag_all_annavg_ens = nan(num_lons,num_lats,length(all_years));

% Calculate ensemble average for each year.
for index=1:length(all_years)
    % Define 3D array, where 3rd dimension is the ensemble member.
    temp_ssh_sep = nan(num_lons,num_lats,num_ens);
    temp_ssh_mar = nan(num_lons,num_lats,num_ens);
%     temp_tilt_mag = nan(num_lons,num_lats,num_ens);
    % Extract "index" time slice from each cell, place it in
    % temp_*(:,:,cell_index)
    for cell_index=1:num_ens %Cycle through each ensemble and copy its data to the temporary arrays defined above.
        temp_ssh_sep(:,:,cell_index) = ssh_all_sep_cell{cell_index,1}(:,:,index);
        temp_ssh_mar(:,:,cell_index) = ssh_all_mar_cell{cell_index,1}(:,:,index);
%         temp_tilt_mag(:,:,cell_index) = tilt_mag_all_cell{cell_index,1}(:,:,index);
    end
    % Calculate the ensemble average by taking the mean of temp_* along the
    % time dimension. Place the resulting average in *_ens(:,:,index).
    ssh_all_annavg_ens_sep(:,:,index) = mean(temp_ssh_sep,3,'omitnan');
    ssh_all_annavg_ens_mar(:,:,index) = mean(temp_ssh_mar,3,'omitnan');
%     tilt_mag_all_annavg_ens(:,:,index) = mean(temp_tilt_mag,3,'omitnan');
end

% SSH climatology for 1950-1970
ssh_all_high_ens_clim_mar = mean(ssh_all_annavg_ens_mar(:,:,high_years_indices),3,'omitnan');
ssh_low_ens_clim_sep = mean(ssh_all_annavg_ens_sep(:,:,low_years_indices),3,'omitnan');
ssh_low_ens_clim_mar = mean(ssh_all_annavg_ens_mar(:,:,low_years_indices),3,'omitnan');
% tilt_mag_all_high_ens_clim = mean(tilt_mag_all_annavg_ens(:,:,steady_years_indices),3,'omitnan');
% tilt_mag_low_ens_clim = mean(tilt_mag_all_annavg_ens(:,:,low_years_indices),3,'omitnan');

if chopped_flag==1 %Only want 1960-2090
    % Now, for each cell element, cut off the first and last 10 years. 
    % Cut off first and last 10 years from each component model
    % Cut off first 10 and last 10 years
    all_years = 1960:2090; % Redefine all_years from 1:131 in previous section so subsetting in this section works.
    all_years = all_years-1949; %First value in all_years is now 11 (neglecting 1-10 to cut them off)

    ssh_all_annavg_ens_sep = ssh_all_annavg_ens_sep(:,:,all_years); 
    ssh_all_annavg_ens_mar = ssh_all_annavg_ens_mar(:,:,all_years);
    % tilt_mag_all_annavg_ens = tilt_mag_all_annavg_ens(:,:,all_years);

    all_years = all_years-10;

    % Extract high and low years
    high_years = 1960:2024; 
    high_years = high_years-1959; %The first index==1 now.
    low_years = 2025:2090; 
    low_years = low_years-1959;
else %Keeping all years and just defining years when speeds increase (high_years) and decrease (low_years; corresponds to our defined "projected" time span). 
    all_years = all_years-1949;
    % Extract high and low years
    high_years = 1950:2024; %OR 2043, if include the flattening out area.
    high_years = high_years-1949; %The first index==1 now.
    low_years = 2025:2100; %OR 2044-2090, if include flattening out area.
    low_years = low_years-1949;    
end

% Extract appropriate data for each year chunk and each ensemble variable:
ssh_low_annavg_ens_sep = ssh_all_annavg_ens_sep(:,:,low_years);
ssh_low_annavg_ens_mar = ssh_all_annavg_ens_mar(:,:,low_years);
ssh_high_annavg_ens_mar = ssh_all_annavg_ens_mar(:,:,high_years);
ssh_high_annavg_ens_sep = ssh_all_annavg_ens_sep(:,:,high_years);
% tilt_mag_low_annavg_ens = tilt_mag_all_annavg_ens(:,:,low_years);
% tilt_mag_high_annavg_ens = tilt_mag_all_annavg_ens(:,:,high_years);

% Calculate trends:
% Transpose low and high years: 
high_years = high_years';
low_years = low_years';
all_years = all_years';

ssh_all_ens_trend_mar = trend_calc(ssh_all_annavg_ens_mar,lat,lon,all_years);
ssh_low_ens_trend_sep = trend_calc(ssh_low_annavg_ens_sep,lat,lon,low_years);
ssh_low_ens_trend_mar = trend_calc(ssh_low_annavg_ens_mar,lat,lon,low_years);
ssh_high_ens_trend_mar = trend_calc(ssh_high_annavg_ens_mar,lat,lon,high_years);
ssh_high_ens_trend_sep = trend_calc(ssh_high_annavg_ens_sep,lat,lon,high_years);
% tilt_mag_all_ens_trend = trend_calc(tilt_mag_all_annavg_ens,lat,lon,all_years);
% tilt_mag_low_ens_trend = trend_calc(tilt_mag_low_annavg_ens,lat,lon,low_years);
% tilt_mag_high_ens_trend = trend_calc(tilt_mag_high_annavg_ens,lat,lon,high_years);

% Already calculated and angle-corrected: u_rhs_geo_*_ens_trend, v_rhs_geo_*_ens_trend
% From section 2, u_rhs_geo_all_high*clim, u_rhs_geo_low*clim, and their v
% counterparts are angle-corrected and already calculated. 

% FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Define function for calculating Arctic-wide speed time series
function out_speed_smoothed = area_avg_speed(three_D_speed_array,weights_of_latitudes,lon_count,lat_count,years,smooth_switch)
    speed_num = zeros(1,length(years)); %Weighted-average speed numerator. Weights based on grid area (uarea, or in the function, "weights_of_latitudes")
    speed_den = zeros(1,length(years)); %Weighted-average speed denominator.
    
    for yr_idx=1:length(years) %Finally, cycle through each year.
        for row_idx=1:lon_count %Then, cycle through longitude.
            for col_idx=1:lat_count %First, cycle through each latitude.
                if ~isnan(three_D_speed_array(row_idx,col_idx,yr_idx)) %Only calculate the numerator and denominators for non-NaN values.
                    speed_num(1,yr_idx) = speed_num(1,yr_idx) + three_D_speed_array(row_idx,col_idx,yr_idx)*weights_of_latitudes(row_idx,col_idx);
                    speed_den(1,yr_idx) = speed_den(1,yr_idx) + weights_of_latitudes(row_idx,col_idx);
                end
            end
        end        
    end
    %Calculate speed_avg array:
    speed_avg = speed_num./speed_den; %Length==number of years.

    if smooth_switch==1 %We want to smoothe the time series within the function.
        %Smooth the data using a 21-year moving window. Use "movmean" command.
        out_speed_smoothed = smoothdata(speed_avg,'movmean',21);
        out_speed_smoothed = out_speed_smoothed(11:141);
    else %We want the unsmoothed output.
        out_speed_smoothed = speed_avg;
    end
end

function final_trend = trend_calc(in_array,latitude_array,longitude_array,year_array)
% Define the number of latitudes and longitudes using latitude_array and
% longitude_array
num_longitudes = size(longitude_array,1);
num_latitudes = size(latitude_array,2);

% Define final_trend as a nan array
final_trend = nan(num_longitudes, num_latitudes);

for row_index=1:num_longitudes %Finally, cycle through each longitude.
    for col_index=1:num_latitudes %First, cycle through each latitude.
        test_years = year_array; %Copy year_array to test_years for the purposes of Matlab regression (i.e., it needs two columns).
        temp_array_gridpoint = in_array(row_index,col_index,:); %Extract a given lat/lon grid point's data.
        temp_array_gridpoint = permute(temp_array_gridpoint,[3 1 2]); %Permute the extracted data so it becomes a column vector.
        temp_nans = isnan(temp_array_gridpoint); %Find which data points are NaNs
        if ~isequal(temp_nans,ones(size(temp_array_gridpoint,1),1)) %Assuming the entirety of temp_array_gridpoint isn't all NaNs, NaN out years in temp_years corresponding to NaNs in temp_array_gridpoint
            for index=1:size(test_years,1) %Cycle through each year.
                if isnan(temp_array_gridpoint(index,1)) %If a given data value is a NaN, the corresponding time must also be NaN'ed. Otherwise, regression won't work.
                    test_years(index,1) = NaN;
                end
            end
            % Now it's time to add the column of 1's so regress actually works:
            test_years(:,2) = test_years(:,1); %Move in_years original data to the newly formed second column
            test_years(:,1) = ones(size(test_years,1),1);
            b = regress(temp_array_gridpoint,test_years); %Calculate the regression. 
            % We're only keeping the slope value, not the y-intercept.
            final_trend(row_index,col_index) = b(2); %slope
        end
        clear temp_array_gridpoint
        clear test_years
        clear temp_nans
    end
end
end
