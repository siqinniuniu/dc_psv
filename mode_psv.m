function [ M, Minv, Q ] = mode_psv( vp, vs, rho, p )
% computes the normalized modal and projection matrices in
% layered elastic isotropic medium for P-SV motion.
% A*d[v]/dt+d[v]/dz=0;  A = M*Q*Minv
%------------------------------------------------------------------
% Usage: [ M, Minv, Q ] = mode_psv( vp, vs, rho, p )
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
%   Q: diagonal matrix of vertical slowness, J = diag([-qp,qp,-qs,qs]);
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
%   [2013-10-20] modified: the inverse matrix is calculated from
%   transpose(M)*I2, as suggested by Chapman<<fundamental of seismic wave propagation>>

%%

% predefined variables
p2 = p^2;
mu = rho*vs^2;
qp = (vp^-2-p2)^(0.5);
qs = (vs^-2-p2)^(0.5);
eta = 2*p2*mu-rho;
trp = 2*p*mu*qp;
trs = 2*p*mu*qs;

% M: matrix of eigenvector ( T: traction vector on z-plane)

%   down-going P      up-going P    down-going SV   up-going SV 
M = [ p,                p,              qs,             qs      ;  % Vr
      qp,              -qp,            -p,              p       ;  % Vz
     -trp,              trp,            eta,           -eta     ;  % Tr
      eta,              eta,            trs,            trs ]   ;  % Tz

% normalize velocity amplitude in M to one     
M = M*diag([vp,vp,vs,vs]); 

% phase matrix
%   down up  down up
Q = [qp; -qp; qs; -qs];

% inverse matrix of M
I2 = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0];
Minv = diag([-qp*vp^2,qp*vp^2,-qs*vs^2,qs*vs^2].^-1)*transpose(M)*I2/rho/2;

end