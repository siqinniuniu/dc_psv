%downward continuate surface plane wavefield by a station on top of the
%ice sheet.

clear all;close all;clc
%% parameters

%time samples
t0 = -10;
t1 = 30;
fs = 20;

%Earth model (ice/bedrock)
vp = [3.87 5.8]; % km/s
vs = [1.95 3.46]; % km/s
rho = [0.917 2.72]; % g/cm^3
% thik = [35 0];  % km
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

%% downward continuate the surface record

nt = length(t);
v0 = zeros(4,nt);
v0(1,:) = vr;
v0(2,:) = -vz;

m1 = dc_psv(...
    nlyr,vp,vs,rho,thik,... % earth model
    nt,v0,fs,...            % time samples of velocity-stress vector
    rayp);                  % ray parameter

%% plot

figure
subplot(3,1,1); plot(t,vr,'k',t,vz,'r',t,rz,'b'); legend('vr','vz','rz')
subplot(3,1,2); plot(t,m1(1,:),t,m1(2,:)); legend('Pd','Pu')
subplot(3,1,3); plot(t,m1(3,:),t,m1(4,:)); legend('Sd','Su')
xlabel('Time (s)')