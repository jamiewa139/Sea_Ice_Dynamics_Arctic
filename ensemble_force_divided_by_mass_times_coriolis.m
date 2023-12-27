% Ensemble VAR/(f*ma) calcluation using daily data. Done separately for
% each ensemble member. VAR options: strair(x/y), strocn(x/y), and
% strint(x/y). Out variables are VAR_d_div_fma. ma==sea ice mass and is
% calculated using rho=917 kg/m3. f==Coriolis parameter, calculated using
% omega=0.73e-4 (1/s). 

% User input values (i.e., in-variable name, out-variable name)
in_dir = 'strinty_d/'; %Options: strinty_d, strinty_d, strinty_d, strinty_d, strinty_d, and strinty_d
num_ens = 10; %Number of ensemble members being analyzed.
root_dir='//wsl$/ubuntu/home/jamiewa/cesm2_lens/';
month_num = 9; %9==September, 3==March
% Corresponding in_var and out_var name info.
in_var = in_dir(1:end-1); %all but end slash is included in the variable name.
out_var = [in_var,'_div_fma'];

% Ensemble directory information (for num_ens)
dir_ens = cell(num_ens,1);
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
% File fronts:
hist_file_front = 'b.e21.BHISTsmbb.f09_g17.LE2-'; %Historical
fut_file_front = 'b.e21.BSSP370smbb.f09_g17.LE2-'; %Future
% File middles:
ice_file_middle = '.cice.h1.'; %Daily sea ice data.
pop_file_middle = '.pop.h.nday1.'; %Daily ocean data.
% Month-based file information:
if month_num==9
    num_days_per_year = 30; %CESM2 uses a 365-day calendar. 
    in_hist_file_end = '.19500101-20141231_arctic_sep.nc'; %input data is daily
    in_fut_file_end = '.20150101-21001231_arctic_sep.nc';
    out_hist_file_end = '.19500101-20141231_arctic_sep_avg.nc'; %output daily is yearly.
    out_fut_file_end = '.20150101-21001231_arctic_sep_avg.nc';
elseif month_num==3
    num_days_per_year = 31; %CESM2 uses a 365-day calendar. 
    in_hist_file_end = '.19500101-20141231_arctic_mar.nc'; %input data is daily
    in_fut_file_end = '.20150101-21001231_arctic_mar.nc';
    out_hist_file_end = '.19500101-20141231_arctic_mar_avg.nc'; %output daily is yearly.
    out_fut_file_end = '.20150101-21001231_arctic_mar_avg.nc';
end

in_thick_dir = 'sithick_d/';
in_thick_var = in_thick_dir(1:end-1); %sithick_d
in_siu_dir = 'siu_d/';
in_siu_var = in_siu_dir(1:end-1); %siu_d

% Define year chunks in each file.
hist_years = 1950:2014;
fut_years = 2015:2100;
num_hist_years = length(hist_years);
num_fut_years = length(fut_years);
% Define physical constants
rho_ice = 917; %kg/m3 (density of ice)
omega = 0.73e-4; %For calculating coriolis
% Import dist2coast data.
dist2coast_data = ncread('//wsl$/ubuntu/home/jamiewa/cmip6_nc_files/dist2coast_cesm2_60deg.nc','dist');

% Import latitude and longitude; calculate the coriolis parameter (f).
lat = ncread([root_dir,in_dir,dir_ens{1,1},hist_file_front,dir_ens{1,1}(1:end-1),ice_file_middle,in_var,in_hist_file_end],'ULAT');
lon = ncread([root_dir,in_dir,dir_ens{1,1},hist_file_front,dir_ens{1,1}(1:end-1),ice_file_middle,in_var,in_hist_file_end],'ULON');
num_lons = size(lon,1);
num_lats = size(lat,2);

% Sometimes CESM2 assigns unrealistically large latitude values where
% variables aren't defined. These lat (and their corresponding longitude)
% values need to be NaN'ed.
for row_idx=1:num_lons
    for col_idx=1:num_lats
        if lat(row_idx,col_idx)>90
            lat(row_idx,col_idx)=NaN;
            lon(row_idx,col_idx)=NaN;
        end
    end
end

% Calculate f (2-dimensional)
f = 2 * omega * sind(lat);

% For each ensemble and for each year, we will import 1) in_var and 2)
% in_thick_var for "num_days_per_year" days. We will calculate ma and then VAR/(f*ma) for
% each day, and then average over all days. This is done separately for
% historical and future files.

% Historical calculations:
for ens_index=1:num_ens %Finally, cycle through each ensemble member.
    % First, predefine output annual average variable array.
    out_hist_annavg = nan(num_lons,num_lats,num_hist_years); %Dimensions: number of longitudes by number of latitudes by (in third dimension) number of historical years.
    disp(['Ensemble Member: ',dir_ens{ens_index,1}])
    day_start_idx = 1; %Counter used to keep track of which day we're starting with during the import process. This is because we're not importing the entire historical file at once.
    day_end_idx = num_days_per_year; %30 or 31, depending on the month.    
    for year_index=1:num_hist_years %1 to 66
        disp(['Importing data for ',num2str(hist_years(year_index)),':', num2str(day_start_idx),' to ', num2str(day_end_idx)])
        start_loc = [1 1 day_start_idx]; %Indices used to define the start and end of time slices to be imported by ncread.
        end_loc = [Inf Inf num_days_per_year];
        % Import the variable of interest for the given time range
        % (start_loc and end_loc)
        temp_in_var(:,:,:) = ncread([root_dir,in_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),...
            ice_file_middle,in_var,in_hist_file_end],in_var,start_loc,end_loc);
        % Import the corresponding sea ice thickness.
        temp_mass_var(:,:,:) = ncread([root_dir,in_thick_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),...
            ice_file_middle,in_thick_var,in_hist_file_end],in_thick_var,start_loc,end_loc);
        % Calculate sea ice mass using rho_ice (==917 kg/m3)
        temp_mass_var(:,:,:) = temp_mass_var(:,:,:)*rho_ice; %Calculate simass from sithick (simass=sithick*rho_ice)
        % Import siu (this is only used to NaN out values in temp_in_var)
        temp_siu_var(:,:,:) = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),...
            ice_file_middle,in_siu_var,in_hist_file_end],in_siu_var,start_loc,end_loc);
        % Now that data has been read in for each ensemble member, it's time to
        % calculate VAR/(f*ma)
        temp_out_var(:,:,:) = temp_in_var(:,:,:)./(f(:,:).*temp_mass_var(:,:,:));
        
        % For all days, NaN out dist2coast cells: 
        for row_idx=1:num_lons %Finally, cycle through each lonigtude.
            for col_idx=1:num_lats %First, cycle through each latitude.
                if dist2coast_data(row_idx,col_idx)<150 || isnan(dist2coast_data(row_idx,col_idx)) %If a given gridpoint is too close to shore or has a NaN dist2coast value, the corresponding variable value is NaN'ed.
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;
                end
            end
        end
        
        % NaN out non-Tandon
        % domain lat/lons (Fig. 1 in Tandon et al., 2018)
        for row_idx=1:num_lons %Finally, cycle through each longitude.
            for col_idx=1:num_lats %first, cycle through each latitude.
                if (lon(row_idx,col_idx)>=103 && lon(row_idx,col_idx)<=236) && lat(row_idx,col_idx)<68 %Latitude has to be greater than 68 degN for longitudes between 103 and 236 degE.
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;
                elseif (lon(row_idx,col_idx)>236 || lon(row_idx,col_idx)<103) && lat(row_idx,col_idx)<79 %Elsewhere, latitude has to be greater than 79 degN.
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;                    
                elseif isnan(lon(row_idx,col_idx)) || isnan(lat(row_idx,col_idx)) %Assuming that NaN lat/lon values slipped through the cracks, this makes sure that corresponding variable gridpoints are NaN'ed out.
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;                    
                end
            end
        end
        
        % Now, for each day, NaN out temp_out_var pixels where temp_siu_var
        % or temp_mass_var==NaN
        for day_index=1:num_days_per_year %Finally, cycle through each year
            for row_idx=1:num_lons %Then, cycle through each longitude.
                for col_idx=1:num_lats %First, cycle through each latitude.
                    % if the corresponing temp_siu_var value is a NaN or
                    % temp_mass_var is a NaN, NaN out temp_out_var.
                    if isnan(temp_siu_var(row_idx,col_idx,day_index)) || isnan(temp_mass_var(row_idx,col_idx,day_index))
                        temp_mass_var(row_idx,col_idx,:)=NaN;
                        temp_siu_var(row_idx,col_idx,:)=NaN;
                        temp_out_var(row_idx,col_idx,:)=NaN;
                    end
                end
            end
        end
            
        clear temp_mass_var % Get rid of the original data
        clear temp_in_var
        % Calculate monthly average of temp_out_var and place 
        out_hist_annavg(:,:,year_index) = mean(temp_out_var,3,'omitnan'); 
        % Clear out temp_out_var (daily data)
        clear temp_out_var
        % Determine the new start day and end day.
        day_start_idx = day_start_idx + num_days_per_year; %Increase start_index by 365 days to move to the next year.
        day_end_idx = day_end_idx + num_days_per_year; %Increase end_index by 365 days.
    end
    % Once we've cycled through all years for a given ensemble member, we
    % write the output in out_hist_annavg(:,:,:) to the
    % appropriate *_MONTH_avg.nc file
    ncwrite([root_dir,in_dir,dir_ens{ens_index,1},hist_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_var,out_hist_file_end],...
        out_var,out_hist_annavg(:,:,:));
end

% Future
for ens_index=1:num_ens %Finally, cycle through each ensemble member.
    % First, predefine output annual average variable array.
    out_fut_annavg = nan(num_lons,num_lats,num_fut_years); %Dimensions: number of longitudes by number of latitudes by (in third dimension) number of historical years.
    disp(['Ensemble Member: ',dir_ens{ens_index,1}])
    day_start_idx = 1; %Define the start index for the first day of the month in the first year of future_years.
    day_end_idx = num_days_per_year; %Depends on the month.    
    for year_index=1:num_fut_years %1 to 86
        disp(['Importing data for ',num2str(fut_years(year_index)),':', num2str(day_start_idx),' to ', num2str(day_end_idx)])
        start_loc = [1 1 day_start_idx]; %Indices used to define the start and end of time slices to be imported by ncread.
        end_loc = [Inf Inf num_days_per_year];
        % Import the variable of interest for a given month.
        temp_in_var(:,:,:) = ncread([root_dir,in_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),...
            ice_file_middle,in_var,in_fut_file_end],in_var,start_loc,end_loc);
        % Import sea ice thickness for a given month.
        temp_mass_var(:,:,:) = ncread([root_dir,in_thick_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),...
            ice_file_middle,in_thick_var,in_fut_file_end],in_thick_var,start_loc,end_loc);
        % Calculate sea ice mass using sea ice thickness and sea ice
        % density values.
        temp_mass_var(:,:,:) = temp_mass_var(:,:,:)*rho_ice; %Calculate simass from sithick (simass=sithick*rho_ice)
        % For future masking purposes only, import siu
        temp_siu_var(:,:,:) = ncread([root_dir,in_siu_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),...
            ice_file_middle,in_siu_var,in_fut_file_end],in_siu_var,start_loc,end_loc);
        % Now that data has been read in for each ensemble member, it's time to
        % calculate VAR/(f*ma)
        temp_out_var(:,:,:) = temp_in_var(:,:,:)./(f(:,:).*temp_mass_var(:,:,:));
        
        % For all days, NaN out dist2coast cells
        for row_idx=1:num_lons %finally, cycle through each longitude.
            for col_idx=1:num_lats %first, cycle through each latitude.
                % If the gridpoint is within 150km of shore or dist2coast
                % is not defined, NaN out sea ice mass, variable, and siu
                % data.
                if dist2coast_data(row_idx,col_idx)<150 || isnan(dist2coast_data(row_idx,col_idx))
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;
                end
            end
        end

        % NaN out appropriate non-Tandon
        % domain lat/lons
        for row_idx=1:num_lons %Finally, cycle through each longitude.
            for col_idx=1:num_lats %First, cycle through each latitude.
                if (lon(row_idx,col_idx)>=103 && lon(row_idx,col_idx)<=236) && lat(row_idx,col_idx)<68 %Between 103 and 236 degE, the minimum latitude is 68 degN
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;
                elseif (lon(row_idx,col_idx)>236 || lon(row_idx,col_idx)<103) && lat(row_idx,col_idx)<79 %At all other longitudes, the minimum latitude is 79 degN
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;                    
                elseif isnan(lon(row_idx,col_idx)) || isnan(lat(row_idx,col_idx)) %If variable information is somehow defined where lat/lon values aren't defined, NaN them out.
                    temp_mass_var(row_idx,col_idx,:)=NaN;
                    temp_siu_var(row_idx,col_idx,:)=NaN;
                    temp_out_var(row_idx,col_idx,:)=NaN;                    
                end
            end
        end
        
        % Now, for each day, NaN out temp_out_var pixels where temp_siu_var
        % or temp_mass_var==NaN
        for day_index=1:num_days_per_year %Finally, cycle through each day
            for row_idx=1:num_lons %Then, cycle through each longitude.
                for col_idx=1:num_lats %First, cycle through each latitude.
                    % If temp_siu_var or temp_mass_var are NaN's at a given
                    % location, NaN out all variable information.
                    if isnan(temp_siu_var(row_idx,col_idx,day_index)) || isnan(temp_mass_var(row_idx,col_idx,day_index))
                        temp_mass_var(row_idx,col_idx,:)=NaN;
                        temp_siu_var(row_idx,col_idx,:)=NaN;
                        temp_out_var(row_idx,col_idx,:)=NaN;
                    end
                end
            end
        end
        
        clear temp_mass_var % Get rid of the original data
        clear temp_in_var
        % Calculate monthly average of temp_out_var and place 
        out_fut_annavg(:,:,year_index) = mean(temp_out_var,3,'omitnan'); 
        % Clear out temp_out_var (daily data)
        clear temp_out_var
        day_start_idx = day_start_idx + num_days_per_year; %Increase start_index by num_days_per_year days to move to the next year.
        day_end_idx = day_end_idx + num_days_per_year; %Increase end_index by num_days_per_year days.
    end
    % Once we've cycled through all years for a given ensemble member, we
    % write the output in out_hist_annavg(:,:,:) to the
    % appropriate *_annavg.nc file
    ncwrite([root_dir,in_dir,dir_ens{ens_index,1},fut_file_front,dir_ens{ens_index,1}(1:end-1),ice_file_middle,in_var,out_fut_file_end],...
        out_var,out_fut_annavg(:,:,:));
end