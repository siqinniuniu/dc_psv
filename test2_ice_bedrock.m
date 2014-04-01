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
m0 = zeros(4*nevt,nt);
m1 = zeros(4*nevt,nt,nz);

rayp = [sacst(:,1).user0];
for ievt = 1:nevt

    v0(1,:) = sacst(ievt,1).data;
    v0(2,:) = -1*sacst(ievt,2).data;
    
    %mode vector in the ice sheet
    idx4_evt = 4*(ievt-1)+(1:4);
    thik = [0 0];
    m0(idx4_evt,:) = dc_psv(...
            1,vp,vs,rho,thik,... % earth model
            nt,v0,fs,...         % time samples of velocity-stress vector
            rayp(ievt));
    
    %mode vector in the bedrock
    for iz = 1:nz
        thik = [z(iz) 0];
        m1(idx4_evt,:,iz) = dc_psv(...
            nlyr,vp,vs,rho,thik,... % earth model
            nt,v0,fs,...            % time samples of velocity-stress vector
            rayp(ievt));            % ray parameter
    end
end

%% calculate Su energy reduction ratio

Esu0 = zeros(nevt,1);
for ievt = 1:nevt
    idx0 = 4*(ievt-1);
    Su = m0(idx0+4,:);
    qs = sqrt(vs(1)^-2-rayp(ievt)^2);
    coef = rho(1)*vs(1)^2*qs;
    Esu0(ievt) = coef*sum(Su.^2);
end

Esu1 = zeros(nevt,nz);
for ievt = 1:nevt
    qs = sqrt(vs(2)^-2-rayp(ievt)^2);
    coef = rho(2)*vs(2)^2*qs;
    
    idx0 = 4*(ievt-1);
    for iz = 1:nz
        Su = m1(idx0+4,:,iz);
        Esu1(ievt,iz) = coef*sum(Su.^2);
    end
end

%Su energy reduction
REDsu = 1-Esu1./repmat(Esu0,1,nz);

%% plot

figure
plot(z,REDsu)
% ylim([-1 1])