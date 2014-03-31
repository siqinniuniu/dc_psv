function [sacst,t] = SACST_synPRF_haskell(...
    nlyr,vp,vs,rho,thik,... % earth model
    t0,t1,fs,...            % time samples
    rayp,a)                 % ray parameter
%compute surface response of plane layered elastic model from 
%plane P-wave incidence (non-evanescent wave)
%sacst(1:3): 1,Vr; 2,Vz; 3,RZ

% clear all;close all;clc
%% parameters

%time samples
% t0 = -10;
% t1 = 30;
% fs = 10;

%Earth model
% vp = [6.4 8.1]; % km/s
% vs = [3.65 4.5]; % km/s
% rho = [2.7 3.3]; % g/cm^3
% thik = [35 0];  % km
% nlyr = 2;

%incident wave
% rayp = 0.06; % s/km

%low-pass filter for making receiver function
% a = 5;
F = @(w,a) exp(-w.^2/4/a^2); 

%% Haskell propagator

%time samples
dt = 1/fs;
t = (t0:dt:2*t1)'; %double the length
nt = length(t);

%frequency samples
dw = 2*pi*fs/nt;
w = dw*(0:nt-1)';
wnyq = pi*fs;
idx = w>(wnyq+dw/4);
w(idx) = w(idx)-2*pi*fs;

%casual exponential taper in time domain to avoid aliasing from frequency sampling
tapercoef = -log(0.01)/2/t1; 
% beta = 0;
w1 = w-1i*tapercoef;

%haskell matrix
H = repmat(eye(4),1,nt); %nt*(4,4) matrix
for ilyr = 1:nlyr-1
    %mode/phase matrix
    [M,Minv,J] = mode_psv(vp(ilyr),vs(ilyr),rho(ilyr),rayp);
    %compute H for each frequency sample
    for iw = 1:nt
        Hilyr = M*diag(exp(1i*w1(iw)*J*thik(ilyr)))*Minv;
        idx4 = 4*(iw-1)+(1:4);
        H(:,idx4) = Hilyr*H(:,idx4);
    end
end

%apply boundary condition in the last layer
[~,Minv] = mode_psv(vp(nlyr),vs(nlyr),rho(nlyr),rayp);
%vr->Pu/Su in the last layer
idx_vr = 1+4*(0:nt-1);
vr2Pu = Minv(2,:)*H(:,idx_vr);
vr2Su = Minv(4,:)*H(:,idx_vr);
%vz->Pu/Su in the last layer
idx_vz = 2+4*(0:nt-1);
vz2Pu = Minv(2,:)*H(:,idx_vz);
vz2Su = Minv(4,:)*H(:,idx_vz);
%time shift the source time function: align direct P to zero time
tp = sqrt(vp.^-2-rayp^2).*thik;
tp = t0+sum(tp(1:nlyr-1));
amp1 = a/pi^0.5*erf(wnyq/2/a);%normalized amplitude in time-domain
S = exp(1i*w1*tp).*F(w1,a)/amp1;
%Pu=Source time function,Su=0 in the last layer
V = zeros(nt,2);
for iw = 1:nt
    V(iw,:) = reshape([vr2Pu(iw),vz2Pu(iw);...
                       vr2Su(iw),vz2Su(iw)]\[S(iw);0],1,[]);
end

%make receiver function
RZ = V(:,1)./V(:,2).*exp(1i*w1*t0).*F(w1,a)/amp1;

%time window cut, calibrate the exponential taper in time domain
idx = t<=t1;
t = t(idx);
%surface velocity
v = fs*ifft(V,'symmetric');
v = v(idx,:).*exp(tapercoef*(t(:,ones(2,1))-t0));
%receiver function
rz = fs*ifft(RZ,'symmetric');
rz = rz(idx).*exp(tapercoef*(t-t0));

%% Output SACST

sacst = SACST_new(1,3);
[sacst.delta] = deal(dt);
[sacst.npts] = deal(length(t));
[sacst.b] = deal(min(t));
[sacst.e] = deal(max(t));
[sacst.user0] = deal(rayp);

sacst(1).kcmpnm = 'R'; sacst(1).data = v(:,1);
sacst(2).kcmpnm = 'Z'; sacst(2).data = -1*v(:,2); %reverse to upward positive
sacst(3).kcmpnm = 'RZ'; sacst(3).data = -1*rz; %Vz:upward positive

end