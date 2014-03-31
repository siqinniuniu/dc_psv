function [ H ] = haskell_dec_psv(...
    nlyr,vp,vs,rho,thik,... % earth model
    nw,w,...                % frequency samples
    rayp)                   % ray parameter
%compute Haskell matrix of plane layered elastic model from 
%plane P-wave incidence (non-evanescent wave)
%size(H) = [4,4*nw];

%% Haskell propagator

%haskell matrix
H = repmat(eye(4),1,nw); %nw*(4,4) matrix
for ilyr = 1:nlyr-1
    %mode/phase matrix
    [M,Minv,J] = mode_psv(vp(ilyr),vs(ilyr),rho(ilyr),rayp);
    %compute H for each frequency sample
    for iw = 1:nw
        Hilyr = M*diag(exp(1i*w(iw)*J*thik(ilyr)))*Minv;
        idx4 = 4*(iw-1)+(1:4);
        H(:,idx4) = Hilyr*H(:,idx4);
    end
end

%decompose wave into (Pd,Pu,Sd,Su) in the last layer
[~,Minv] = mode_psv(vp(nlyr),vs(nlyr),rho(nlyr),rayp);
for iw = 1:nw
    idx4 = 4*(iw-1)+(1:4);
    H(:,idx4) = Minv*H(:,idx4);
end

end