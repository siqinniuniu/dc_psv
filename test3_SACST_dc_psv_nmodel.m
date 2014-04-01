%downward continuate surface plane wavefield by a station on top of the
%ice sheet.

clear all;close all;clc
%% parameters

%sac
sacdir = 'sac/';
saclst = [sacdir,'evt.lst'];
cmpnm = {'.BHR','.BHZ'};

%time samples
t0 = -5;
t1 = 10;
fs = 20;

%only vary thickness of ice sheet
z = 1:0.1:5;
nmod = length(z);

%Earth model (ice/bedrock)
vp0 = [3.87; 5.8]; % km/s
vs0 = [1.95; 3.46]; % km/s
rho0 = [0.917; 2.72]; % g/cm^3
nlyr = 2;

%model space
vp = repmat(vp0,1,nmod);
vs = repmat(vs0,1,nmod);
rho = repmat(rho0,1,nmod);
thik = zeros(nlyr,nmod); thik(1,:) = z;

% sachd
sachd_rayp = 'user0';

%% downward continuation

[REDsu,m1,Esu1,m0,Esu0] = SACST_dc_psv_nmodel(...
    sacdir,saclst,cmpnm,...        % sac data
    t0,t1,fs,...                   % time window for downward continuation
    nlyr,nmod,vp,vs,rho,thik,...   % earth models 
    sachd_rayp);                    % head field of ray parameter

%% plot

figure
plot(z,REDsu)
% ylim([-1 1])