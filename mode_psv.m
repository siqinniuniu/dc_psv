function [ M, Minv, J ] = mode_psv( vp, vs, rho, p )
% computes the normalized modal and projection matrixes
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
      qp,              -qp,            -p,              p			;  % Vz
     -2*mu*p*qp,        2*mu*p*qp,     -mu*qss,         mu*qss		;  % Tr
     -mu*qss,          -mu*qss,         2*mu*p*qs,      2*mu*p*qs	]; % Tz

% modal projection matrix
%        Vr               Vz                 Tr          Tz
Minv = [ p*mu,            mu*qss/qp/2,      -p/qp/2,    -0.5       ;  % down-going P
         p*mu,           -mu*qss/qp/2,       p/qp/2,    -0.5       ;  % up-going P
         mu*qss/qs/2,    -p*mu,             -0.5,		 p/qs/2    ;  % down-going SV
         mu*qss/qs/2,     p*mu,              0.5,		 p/qs/2    ]/rho;  % up-going SV

% normalize velocity amplitude in M to one
M = M*diag([vp,vp,vs,vs]);
Minv = diag([vp^-1,vp^-1,vs^-1,vs^-1])*Minv;

% phase matrix
%   down  up  down up
J = [-qp, qp, -qs, qs];

end