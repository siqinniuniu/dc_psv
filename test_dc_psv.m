%make synthetic plane P-wave response in plane layered elastic model and
%downward continuate

clear all;close all;clc
%% parameters

%time samples
t0 = -10;
t1 = 30;
fs = 20;

%Earth model
vp = [6.4 8.1]; % km/s
vs = [3.65 4.5]; % km/s
rho = [2.7 3.3]; % g/cm^3
thik = [35 0];  % km
nlyr = 2;

%incident wave
rayp = 0.06; % s/km

%low-pass filter for making receiver function
a = 5;

%% SACST_synPRF_haskell

[sacst,t] = SACST_synPRF_haskell(nlyr,vp,vs,rho,thik,t0,t1,fs,rayp,a);

vr = sacst(1).data;
vz = sacst(2).data;
rz = sacst(3).data;

% figure
% plot(t,vr,'k',t,vz,'r',t,rz,'b')

%% downward continuate the surface record

nt = length(t);
v0 = zeros(4,nt);
v0(1,:) = vr;
v0(2,:) = -vz;

m1 = dc_psv(...
    nlyr,vp,vs,rho,thik,... % earth model
    nt,v0,fs,...            % time samples of velocity-stress vector
    rayp);                  % ray parameter