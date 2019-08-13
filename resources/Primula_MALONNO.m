clear all, clc, close all, fs = 12; lw = 2; fn = 'Times New Roman';
tic
% <PRIMULA> 

% Cislaghi et al. (2017)
% Including root reinforcement variability in ...
% Earth Surface Processes and Landforms, 42(12), 1789�1806.

% Cislaghi et al. (2018)
% A probabilistic multidimensional approach to quantify ...
% Geomorphology, 306, 108�127.

% Do you want to display the figures? yes -> fig = 1;
fig = 0;

% Number of Monte Carlo simulations
N = 1000; % Cislaghi et al., 2018
% N = 2000; % Hammond et al., 1992

% Set the path of the folder
Path = '~/PRIMULA/PRIMULA/';
% Set the path of the folder in which the maps of catchment are located
Folder = 'mappe/';
Output = 'Output_';

% Load the raster 
% -> DEM-derived maps obtained using the functions of GRASS into QGIS

% Digital elevation model (DEM)
[DEM,R_DEM] = geotiffread([Path,Folder,'dem/dem_MALONNO_v2.tif']);

% Slope (rad)
[SLOPE,Rinfo] = geotiffread([Path,Folder,'dem/slope_MALONNO.tif']);
SLOPE(SLOPE == -99999) = NaN;
geoinfo = geotiffinfo([Path,Folder,'dem/slope_MALONNO.tif']);

% check!
% figure, mapshow(SLOPE,Rinfo,'DisplayType','surface'); 

% Topographic wetness index (-)
[TWI,R] = geotiffread([Path,Folder,'dem/twi_MALONNO.tif']);

% -> DEM-derived maps of distances from the nearest stream
% used for the assessment of the probability that an amount of sediment or 
% large wood reachs the nearest stream (Cislaghi et al. 2018)

% Maximum flow path (m)
[MFP,R] = geotiffread([Path,Folder,'dem/maximum_flow_path_MALONNO.tif']);

% Overflow path (m)
[OFP,R] = geotiffread([Path,Folder,'dem/overflow_distance_MALONNO.tif']);

% Sub-model PROBSLOPE
% -> filter of the slope map 
% it is possible that the presence of retaining wall or exposed rock faces 
% affects the calculation of the slope (too steep)
 
% Load the map of PROBSLOPE (output of the sub-model)
[PROBSLOPE,R] = geotiffread([Path,Folder,'dem/PROBSLOPE_MALONNO.tif']);
PROBSLOPE(isnan(SLOPE)) = NaN;
PROBSLOPE(PROBSLOPE<0) = 1.1028656e-06;

% check!!
% figure, mapshow(rad2deg(PROBSLOPE),Rinfo,'DisplayType','surface'); 
% figure, mapshow(rad2deg(SLOPE),Rinfo,'DisplayType','surface'); 
% figure, mapshow(rad2deg(SLOPE-PROBSLOPE),Rinfo,'DisplayType','surface'); 

% Map of soils (scale 1:200,000)
% -> http://www.geoportale.regione.lombardia.it/download-dati
% -> freely available 
% (Carta Litologica Regione Lombardia)

% From the Map of soils -> rasterize in function of soil texture
[SOILS,R] = geotiffread([Path,Folder,'dem/soils_MALONNO_v2.tif']);

% From the Map of soils -> rasterize in function of maximum soil depth 
% (4 classes)
[SOILS_z,R] = geotiffread([Path,Folder,'dem/soils_MALONNO.tif']);
S = shaperead([Path,Folder,'sito_studio/Pedologia_25k_MALONNO.shp']);

% Extraction of input parameter values from the probability distribution
% function

% Internal friction angle, Unit weight of soil, soil cohesion, depth...
Soil            = struct();
Soil.Id         = unique([S.COD_UTS1]);
for i = 1:length([Soil.Id])
    % soil depth (m) <- triangular distribution (Hammond et al. 1992)
    ind = find([S.COD_UTS1] == Soil.Id(i),1,'first');
    Soil.maxZ(i) = S(ind).PROF_UTILE ./ 100;
    Soil.Z(i,1:N) = icdf(makedist('Triangular','a',2./3.*Soil.maxZ(i),'b',3./4.*Soil.maxZ(i),'c',Soil.maxZ(i)), rand(N,1));  
end

% check!
% figure, mapshow(SOILS,R,'DisplayType','surface'); caxis([0 2])

% unit weight of the water (kN/m3)
gamma_w = 9.81;

% internal friction angle <- uniform distribution
% http://www.geotechdata.info/parameter/angle-of-friction.html
% set the lower and upper bounds
phi1 = icdf(makedist('Uniform','lower',30,'upper',40), rand(N,1)); % deg
phi2 = icdf(makedist('Uniform','lower',35,'upper',40), rand(N,1)); % deg

% internal friction angle <- uniform distribution
% http://www.geotechdata.info/parameter/angle-of-friction.html
% set the lower and upper bounds
gamma1 = icdf(makedist('Uniform','lower',17,'upper',19), rand(N,1)); % kN/m3
gamma2 = icdf(makedist('Uniform','lower',15,'upper',18), rand(N,1)); % kN/m3

% Cumulative distribution function of soil cohesion (negligible)
cs1 = icdf(makedist('Uniform','lower',0,'upper',1), rand(N,1)); %kPa
cs2 = icdf(makedist('Uniform','lower',0,'upper',1), rand(N,1)); %kPa

% soil permeability coefficient  <- uniform distribution
% http://www.geotechdata.info/parameter/permeability.html
% set the lower and upper bounds
ks1 = icdf(makedist('Uniform','lower',0.5,'upper',100), rand(N,1)); % m/day
ks2 = icdf(makedist('Uniform','lower',0.5,'upper',100), rand(N,1)); % m/day

% Size of landslides <- normal distribution
% Landslide size from the literature (Milledge et al. 2014)
area = 10.^ icdf(makedist('Normal','mu',2.017,'sigma',0.176466^0.5), rand(N,1)); % m2
l2w = 10.^ icdf(makedist('Normal','mu',0.1528,'sigma',0.037396^0.5), rand(N,1)); % dimensionless

% Extraction of width and length of potential landslide
widthCELL = sqrt(area./l2w);        % m
lengthCELL = widthCELL .*  l2w;     % m

% Map of land-use (DUSAF 5.0)
% -> http://www.geoportale.regione.lombardia.it/download-dati
% -> freely available 

% Land use -> rasterize (each code describes a land use)
[DUSAF,R] = geotiffread([Path,Folder,'dem/dusaf_MALONNO.tif']);

% Empirical data of root reinforcement
% For each different forest area, PRIMULA requests different probability
% distribution function that depends of species, density, DBH, etc.

% These data will be available for several forest categories before the end
% of the project
Cr_forest = readtable([Path,Folder,'rinforzo_radicale/RootReinforcement.csv']); % kPa

% Empirical distribution for different forest population density
n = randi(length(Cr_forest.Pa200),N,1);
Crl_Fs200 = Cr_forest.Fs200(n);
Crl_Fs800 = Cr_forest.Fs800(n);
Crl_Pa200 = Cr_forest.Pa200(n);
Crl_Pa400 = Cr_forest.Pa400(n);
Crl_Mf300 = Cr_forest.MF300(n);
Crl_Mf600 = Cr_forest.MF600(n);
Crl_Cs150 = Cr_forest.Cs150(n);

% Cumulative distribution function of Root reinforcement for "Grassland"
% from the literature for depth less than 0.4 m
Cr_grassland = icdf(makedist('Uniform','lower',5,'upper',7.5), rand(N,1));

% Cumulative distribution function of Root reinforcement for "Shrubland"
% from the literature for depth less than 0.4 m
Cr_shrubland = icdf(makedist('Uniform','lower',0,'upper',15), rand(N,1));

% PRIMULA -> geotechnical part
% Initialization of the raster for each parameter

% Cislaghi et al. 2017
% Crb = double(SOILS) * 0;                % Basal root reinforcement (kPa)
% Crl = double(SOILS) * 0;                % Lateral root reinforcement (kPa)
% phi = double(SOILS) * 0;                % Soil friction angle (rad)
% gamma_s = double(SOILS) * 0;            % Specific weight values (kN/m3)
z = double(SOILS) * 0;                  % Soil depth (m)
w = double(SOILS) * 0;                  % Width of cell (m)
l = double(SOILS) * 0;                  % Length of cell (m)
% ks = double(SOILS) * 0;                 % Permeability of soil (m/day)
% m = double(PEDOLOGIA) * 0;            % Soil wetness (-)
% FS = double(SOILS) * 0;                 % Factor of safety (-)
Pr_of_failure = double(SOILS) * 0;      % Probability of failure (-)

% Cislaghi et al. 2018
% Pr. that potential landslide runout reaches the infrastructure (-)
Pr_hazard_runout = double(SOILS) * 0;   

% Probability that potential landslide runout reaches the channel network (-)
Pr_hazard_channel = double(SOILS) * 0;  

% TOPMODEL + LSPP (http://idro.arpalombardia.it/manual/dati_link.html)
% Work in progress...
% Several hydrological models are being tested...
rainfall = double(SOILS) * 0;           % Design rainfall (m/day)
rainfall = rainfall + 0.100;            % m/day
% Tr5   = 0.100 m/day
% Tr50  = 0.175 m/day
% Tr100 = 0.200 m/day
% Tr200 = 0.230 m/day

% waitbar
f = waitbar(0,'Please wait...');

for i = 1:N

     waitbar(i/N,f,sprintf('%d',i))
    
    % Initialization of the raster for each parameter
    % for each simulation, each raster has to be initializated
    FS = double(SOILS) * 0;         % Factor of safety (-)
    phi = double(SOILS) * 0;        % Soil friction angle (rad)
    gamma_s = double(SOILS) * 0;    % Specific weight values (kN/m3)
    ks = double(SOILS) * 0;         % Permeability of soil (m/day)
    Crl = double(SOILS) * 0;        % Lateral root reinforcement (kPa)
    Crb = double(SOILS) * 0;        % Basal root reinforcement (kPa)

    % close all figures!
    close all
    
    % Set the parameters values for the input raster
    phi(SOILS == 1)     = degtorad(phi1(i));
    phi(SOILS == 2)     = degtorad(phi2(i));
    ks(SOILS == 1)      = degtorad(ks1(i));
    ks(SOILS == 2)      = degtorad(ks2(i));
    
    gamma_s = gamma_s + gamma1(i);
    
    for k = 1:length([Soil.Id])
        z(SOILS_z == Soil.Id(k)) = Soil.Z(k,i);
        % figure, mapshow(z,R,'DisplayType','surface');       
    end
        
    % check!
    % figure, mapshow(SOILS,R,'DisplayType','surface');       
    % figure, mapshow(phi,R,'DisplayType','surface');       
    % figure, mapshow(z,R,'DisplayType','surface');  
    
    % Hydrological model (TOPMODEL)
    % m <- saturation degree (m/m)
    m = TOPMODEL_v3(rainfall,ks,z, TWI,SLOPE);
    
    % check!
    % figure, mapshow(m,R,'DisplayType','surface'); caxis([0 1])

    % The watershed is completly covered by pasture and forest
    % 3211, 3212, ... are codes of the DUSAF 5.0
    
    % Lateral root reinforcement (kPa)
    Crl(DUSAF == 3211)  = Cr_grassland(i);
    Crl(DUSAF == 3212)  = Cr_grassland(i);
    Crl(DUSAF == 3221)  = Cr_grassland(i);
    Crl(DUSAF == 332)   = Cr_shrubland(i);
    Crl(DUSAF == 333)   = Cr_shrubland(i);
    Crl(DUSAF == 3121)  = Crl_Pa400(i);
    Crl(DUSAF == 3122)  = Crl_Pa200(i);
    Crl(DUSAF == 31111) = Crl_Fs800(i);
    Crl(DUSAF == 31121) = Crl_Fs200(i);
    Crl(DUSAF == 3114)  = Crl_Cs150(i);
    Crl(DUSAF == 222)   = Crl_Cs150(i);
    Crl(DUSAF == 31311) = Crl_Mf600(i);
    Crl(DUSAF == 31321) = Crl_Mf300(i);
    

    % Basal root reinforcement (kPa)
    Crb(z >= 0.5)       = 0;
    Crb(z < 0.5)        = Crl(z < 0.5);
       
    % check!
    % figure, mapshow(Crl,R,'DisplayType','surface');
    % figure, mapshow(Crb,R,'DisplayType','surface');
    % figure, mapshow(z,R,'DisplayType','surface');

    w(:) = widthCELL(i);
    l(:) = lengthCELL(i);
    
    % MD-STAB (Milledge et al., 2014) -> application
    FS = MDSTAB_v2(PROBSLOPE,SLOPE,phi,m,gamma_s,w,z,l,Crl,Crb);
    
    % Updating of Pr_of_failure
    Pr_of_failure = Pr_of_failure + double(FS<1);
    
    % check!
    % figure, mapshow(Pr_of_failure,R,'DisplayType','surface');
    
    % Estimation of landslide runout
    Runout = double(FS<1);
    
    % Empirical formulation (Cislaghi et al. 2018)
    Runout = Runout .* 22.82 .* (z .* w .* l) .^ 0.3616;
    
    % Updating of Pr_hazard_runout and Pr_hazard_channel
    Pr_hazard_runout = Pr_hazard_runout + double(Runout>MFP);
    Pr_hazard_channel = Pr_hazard_channel + double(Runout>OFP);
    
end

toc

delete(f)

% check! Probability distribution of factor of safety

% Pr_of_failure is a new variable representing the probability of failure
Pr_of_failure = Pr_of_failure ./ N;
Pr_hazard_runout = Pr_hazard_runout ./ N;
Pr_hazard_channel = Pr_hazard_channel ./ N;
Pr_of_failure(isnan(PROBSLOPE)) = NaN;

% check!
figure, mapshow(Pr_of_failure,R,'DisplayType','surface');
map = [ [0:0.1:1, repmat(1,1,10)]' [repmat(1,1,10) 1:-0.1:0]' [repmat(0,1,21)]' ];
colormap(map)
% figure, mapshow(Pr_hazard_runout,R,'DisplayType','surface');
% figure, mapshow(Pr_hazard_channel,R,'DisplayType','surface');

% check!
% figure, mapshow(Pr_of_failure,R,'DisplayType','surface'); caxis([0 1]);
% color(5,:) = [0 0.5 0];
% color(4,:) = [0 0.92 0];
% color(3,:) = [1 1 0];
% color(2,:) = [1 0.6 0];
% color(1,:) = [0.92 0 0];
% colormap(flipud(color));
% cbh = colorbar;
% hl = ylabel(cbh,'Probability [FS<1] (-)');
% hl.Rotation = 270; hl.Position = [3 0.5 0];
% hl.FontSize = fs;
% print('Figura_finale.png','-dpng','-r300')

% save the result as raster
PrF = Pr_of_failure;
PrF(isnan(SLOPE)) = -999;
t = Tiff([Path,Output,'.tif'], 'w');
Prova_filename = [Path,Output,'.tif'];
geotiffwrite(Prova_filename, PrF, R, 'GeoKeyDirectoryTag', geoinfo.GeoTIFFTags.GeoKeyDirectoryTag)

% TOPMODEL

function [ W ] = TOPMODEL_v3( rainfall, ks, z, TWI, alpha)

% TOPMODEL 
% Hydrological model -> topographic approach

rainfall = rainfall.*0.9; % effects of vegetation interception/evapotranspiration

W = ( rainfall ./ (ks .* z .* cos(alpha)) ) .* TWI;
W(W>1) = 1;

% check!
% figure, mapshow(TWI,R,'DisplayType','surface'); 
% figure, mapshow(W,R,'DisplayType','surface'); caxis([0 1])
end

