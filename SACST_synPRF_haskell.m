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
    [M,Minv,J] = modePSV(vp(ilyr),vs(ilyr),rho(ilyr),rayp);
    %compute H for each frequency sample
    for iw = 1:nt
        Hilyr = M*diag(exp(1i*w1(iw)*J*thik(ilyr)))*Minv;
        idx4 = 4*(iw-1)+(1:4);
        H(:,idx4) = Hilyr*H(:,idx4);
    end
end

%apply boundary condition in the last layer
[~,Minv] = modePSV(vp(nlyr),vs(nlyr),rho(nlyr),rayp);
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

%% nested function
    function [ M, Minv, J ] = modePSV( vp, vs, rho, p )
        % ModePSV computes the normalized modal and projection matrixes
        %   A = M*J*Minv
        %------------------------------------------------------------------
        % Usage: [ M, Minv, J ] = ModeMatrixPSV( vp, vs, rho, p )
        %------------------------------------------------------------------
        % Inputs:
        %   vp: P wavespeed (1 by 1 scalar)
        %   vs: S wavespeed (1 by 1 scalar)
        %   rho: density (1 by 1 scalar)
        %   p: horizontal slowness (1 by 1 scalar)
        %------------------------------------------------------------------
        % Outputs:
        %   M: modal matrix (4 by 4 matrix),
        %       M(:,1|3): down-going P|SV, M(:,2|4): up-going P|SV;
        %   Minv: modal projection matrix (4 by 4 matrix);
        %       Minv(1|3,:): down-going P|SV, Minv(2|4,:): up-going P|SV.
        %   J: diagonal matrix of vertical slowness, J = diag([-qp,qp,-qs,qs]);
        %------------------------------------------------------------------
        % Notice:
        % 1) for near critical rayp (qs or qp = 0), 1/0 will occur in Minv.
        %	However this will be circumvented by combining the phase delay term
        %	in the propagation matrix which will lead to the terms like: sin(qs*tau)/qs
        %	and no singularity occurs at qs=0;
        % 2) vertical direction points downward to the earth's center;
        %------------------------------------------------------------------
        %Appendix
        % Polarization of down- and up-going waves
        %
        %                     Down-going |  Up-going
        %  (tangential) T.-------------------------------> R(radial)
        %                |      SV       |       P
        %                |     /         |      /
        %                | SH *          |  SH *
        %                |     \         |      \
        %                |       P       |       SV
        %                V
        %                Z(down)
        %------------------------------------------------------------------
        % History
        %   [2012-05-27] created
        %   [2012-06-09] modified: change the polarization
        %   [2012-11-17] modified: change the polarization
        
        %%
        
        % predefined variables
        p2 = p^2;
        vs2 = vs^2;
        mu = rho*vs2;
        qp = (vp^-2-p2)^(0.5);
        qs = (1/vs2-p2)^(0.5);
        qss = 1/vs2-2*p2;
        
        %     down-going P      up-going P      down-going SV   up-going SV
        M = [ p,                p,              qs,             qs			;  % Vr
            qp,               -qp,            -p,             p			;  % Vz
            -2*mu*p*qp,       2*mu*p*qp,      -mu*qss,        mu*qss		;  % Tr
            -mu*qss,          -mu*qss,        2*mu*p*qs,      2*mu*p*qs	]; % Tz
        
        % modal projection matrix
        %        Vr               Vz                 Tr          Tz
        Minv = [ p*mu,            mu*qss/qp/2,       -p/qp/2,	 -0.5       ;  % down-going P
            p*mu,            -mu*qss/qp/2,      p/qp/2,	 -0.5       ;  % up-going P
            mu*qss/qs/2,     -p*mu,             -0.5,		 p/qs/2     ;  % down-going SV
            mu*qss/qs/2,     p*mu,               0.5,		 p/qs/2     ]/rho;  % up-going SV
        
        % normalize velocity amplitude in M to one
        M = M*diag([vp,vp,vs,vs]);
        Minv = diag([vp^-1,vp^-1,vs^-1,vs^-1])*Minv;
        
        % phase matrix
        %   down  up  down up
        J = [-qp, qp, -qs, qs];
        
    end


end