% MULTIDIMENSIONAL SHALLOW LANDSLIDE MODEL
%% MD-STAB (Milledge et al., 2014)
% Written by Cislaghi (2015)
% MD-STAB satisfies horizontal and vertical force equilibrium ignoring
% moment equilibrium

% teta      <- slope inclination [rad]
% delta     <- angle from horizontal [rad]
% phi       <- effective friction angle of soil [rad]
% m         <- saturation ration (= hw / z) [-]
% gamma_s   <- unit weigth of the soil ( = rho_s * g)
% rho_s     <- constant bulk density of soil [Kg/m3]
% w         <- width of block [m]
% z         <- vertical depth [m]
% l         <- length of block [m]
% Crl       <- lateral cohesion [Kpa] (depth averaged cohesion) 
% Crb       <- basal cohesion [Kpa]

function [FS,Frd,Frb,Frl,Fdu,Fdc,Kp,Ka] = MDSTAB_v2(teta,delta,phi,m,gamma_s,w,z,l,Crl,Crb)
gamma_w = 9.81; % Unit weigth of water [KN/m3]

%% Central Block Driving Force [KN]
% Equation 3
Fdc = gamma_s .* z .* sin(teta) .* cos(teta) .* w .* l ; 

%% Block Cross-Slope Boundaries [KN]
K0 = 1 - sin(phi);
% Equation 7
Frl = 0.5 .* K0 .* (gamma_s - gamma_w .* m.^2) .* l .* z.^2 .* cos(teta) .* tan(phi) + Crl .* l .* z .* cos(teta);

%% Block Upslope (Fdu) 
% Rankine solution for cohesive soils (Mazindrani and Ganjali, 1997)
% K is the part of equation stated under the square root
K = 4 .* cos(teta).^2 .* (cos(teta).^2 - cos(phi).^2) + (4 .* ((Crl./(gamma_s.*z)).^2) .*cos(phi).^2) + (8 .* (Crl./(gamma_s.*z)) .* cos(teta).^2 .*sin(phi) .*cos(phi));
% Assumption: when slope is greater than teta and K is less than 0
% K assumes the value equal to zero
K(K<0) = 0; 
% Active earth pressure coefficient
% Equation 8
Kp = (1 ./ (cos(phi).^2)) .* (2 .*cos(teta).^2 + 2 .* (Crl./(gamma_s.*z)) .*cos(phi) .*sin(phi) + sqrt(K)) - 1;
% Passive earth pressure coefficient
% Equation 8
Ka = (1 ./ (cos(phi).^2)) .* (2 .*cos(teta).^2 + 2. * (Crl./(gamma_s.*z)) .*cos(phi) .*sin(phi) - sqrt(K)) - 1;
% Net driving force on the upslope margin (slope-parallel component)
% For soils with a strong cohesive component Fdu is negative since the
% resisting forces due to cohesion exceed the driving force of the upslope
% wedge
Fdu = 0.5 .* Ka .* z.^2 .*(gamma_s - gamma_w .* m.^2) .* w .* cos(delta-teta);
% Net driving force on the upslope margin (normal component)
Fnu = 0.5 .* Ka .* z.^2 .*(gamma_s - gamma_w .* m.^2) .* w .* sin(delta-teta);

%% Downslope Boundaries (Frd) 
% Passive force on the downslope margin (slope-parallel component)
Frd = 0.5 .* Kp .* z.^2 .*(gamma_s - gamma_w .* m.^2) .* w .* cos(delta-teta);
% Frd is negligible
Frd = 0;
% Passive force on the downslope margin (normal component)
Fnd = 0.5 .* Kp .* z.^2 .*(gamma_s - gamma_w .* m.^2) .* w .* sin(delta-teta);

%% Basal Resistence Force (Frb)
% Fnt is the effective normale stress on the failure surface integrated
% over the area = Fnc + Fnu - Fnd
Fnc = (gamma_s - gamma_w .*m) .* z .* cos(teta).^2 .* w .* l;
% Equation 14
Fnt = Fnc + Fnu - Fnd;
% Equation 15
Frb = Crb .* w .* l + Fnt .* tan(phi);

%% Factor of Safety
% Equation 16
FS = (Frb + 2. * Frl + Frd  -Fdu) ./ Fdc;
end

