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

%Earth model (ice/bedrock)
vp = [3.87 5.8]; % km/s
vs = [1.95 3.46]; % km/s
rho = [0.917 2.72]; % g/cm^3
% thik = [35 0];  % km
nlyr = 2;

%thickness of ice sheet
z = 1:0.1:5;

%% load SAC, time window cut

sacst = SACST_fread('list',saclst,'prefix',sacdir,'suffix',cmpnm);

dt = 1/fs;
t = t0:dt:t1;
sacst = SACST_interp(sacst,t,'ref','0','rmnan',1);

%% downward continuate the surface wavefield

nz = length(z);
nevt = size(sacst,1);
nt = length(t);
v0 = zeros(4,nt);
m0 = zeros(4*nevt,nt,nz);
m1 = m0;

rayp = [sacst(:,1).user0];
for ievt = 1:nevt
    
    %state vector on the surface
    v0(1,:) = sacst(ievt,1).data;
    v0(2,:) = -1*sacst(ievt,2).data;
    
    %mode vector in the ice/bedrock
    idx4_evt = 4*(ievt-1)+(1:4);
    for imod = 1:nmod
        [mm1,mm0] = dc_psv(...
            nlyr,vp(:,imod),vs(:,imod),rho(:,imod),thik(:,imod),... % earth model
            nt,v0,fs,...            % time samples of velocity-stress vector
            rayp(ievt));            % ray parameter
        m1(idx4_evt,:,iz) = mm1;
        m0(idx4_evt,:,iz) = mm0;
    end
end

%% calculate Su energy reduction ratio

Esu0 = zeros(nevt,nz);
Esu1 = Esu0;
for ievt = 1:nevt
    
    qs0 = sqrt(vs(1)^-2-rayp(ievt)^2);
    coef0 = rho(1)*vs(1)^2*qs;
    
    qs1 = sqrt(vs(nlyr)^-2-rayp(ievt)^2);
    coef1 = rho(nlyr)*vs(nlyr)^2*qs;
    
    idx0 = 4*(ievt-1);
    for iz = 1:nz
        Su1 = m1(idx0+4,:,iz);
        Esu1(ievt,iz) = coef*sum(Su1.^2);
        Su1 = m1(idx0+4,:,iz);
        Esu1(ievt,iz) = coef*sum(Su1.^2);
    end
end

%Su energy reduction
REDsu = 1-Esu1./Esu0;

%% plot

figure
plot(z,REDsu)
% ylim([-1 1])