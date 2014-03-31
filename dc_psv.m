function [ m1 ] = dc_psv(...
    nlyr,vp,vs,rho,thik,... % earth model
    nt,v0,fs,...            % time samples of velocity-stress vector
    rayp)                   % ray parameter
%downward continuate surface plane wavefield through plane layered elastic 
%model (non-evanescent wave), and decompose in the last layer
%size(v0) = size(v0) = [4,nt];

%% parameters

%time shift of S wave during downward continuation
ts = sqrt(vs.^-2-rayp^2).*thik;
ts = sum(ts(1:nlyr-1));

%pad zeros to avoid time warp
npad = ceil(ts*fs);
nw = nt+npad;

%frequency samples
dw = 2*pi*fs/nw;
w = dw*(0:nw-1);
wnyq = pi*fs;
idx = w>(wnyq+dw/4);
w(idx) = w(idx)-2*pi*fs;

%haskell matrix
Hdec  = haskell_dec_psv(...
    nlyr,vp,vs,rho,thik,... % earth model
    nw,w,...                % frequency samples
    rayp);                  % ray parameter

%pad zeros and fft
v0 = [v0,zeros(4,npad)];
Fv0 = fft(v0,[],2);

%apply haskell matrix to surface wavefield
Fm1 = zeros(size(v0));
for iw = 1:nw
    idx4 = 4*(iw-1)+(1:4);
    Fm1(:,iw) = Hdec(:,idx4)*Fv0(:,iw);
end

%apply time shift, ifft and cut
m1 = ifft(repmat(exp(1i*w*ts),4,1).*Fm1,[],2,'symmetric');
m1 = m1(:,1:nt);

end