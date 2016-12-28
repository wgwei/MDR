function [p2lvlSimpd, p2lvlTheo, L] = double_diffr(spos, rpos, N1p, N2p, fr, betas, betar)
% calculate the attenuation of 2.5D thick barrier

 [Ys, Yr, Ygr, Ysm, L]  = getDoubleDiffrPara(spos, rpos, N1p, N2p, fr, betas, betar);
 pSquare = ((0.37./(Ygr+0.37)).^2).*((0.37./(Ysm+0.37)).^2);
 ptheo2 = (f(Ygr)^2+g(Ygr)^2)*(f(Ysm)^2+g(Ysm)^2);
 p2lvlSimpd =  10.0*log10(pSquare); 
 p2lvlTheo = 10*log10(ptheo2);
end

function vf = f(x)
    [c, s] = fcs(x);
    vf  = (0.5-s)*cos(0.5*pi*x^2)-(0.5-c)*sin(0.5*pi*x^2);
end

function vg = g(x)
    [c, s] = fcs(x);
    vg = (0.5-c)*cos(0.5*pi*x^2)+(0.5-s)*sin(0.5*pi*x^2);
end 

 function [Ys, Yr, Ygr, Ysm, L]    = getDoubleDiffrPara(spos, rpos, N1p, N2p, fr, betas, betar)
% calculate the parameters used for double diffraction
% this script is only valide for 2.5D thick barrier
% the parameters are based on A.D.Pierce: Diffraction of sound around
% corners and over wide barriers, 1976, JASA
% [Ys, Yr, Ygr, Ysm, L] = getDoubleDiffrPara(spos, rpos, N1p, N2p, fr)
%   N1p = the upper left point of the barrier
%   N2p = the upper right point of the barrier


rs = dist2D(spos, N1p);
angles = asin(abs(N1p(1)-spos(1))/rs);
thetas = betas+angles-1.5*pi;

rr = dist2D(rpos, N2p);
angler = asin(abs(rpos(1)-N2p(1))/rr);
thetar = betar+angler-1.5*pi;

w = abs(N2p(1)-N1p(1));
L = rs+w+rr;

nus = pi/betas;
nur = pi/betar;
Mvs = (cos(nus*pi)-cos(nus*(betas-thetas)))/(nus*sin(nus*pi));
Mvr = (cos(nur*pi)-cos(nur*(betar-thetar)))/(nur*sin(nur*pi));

Ys = zeros(length(fr), 1);
Yr = zeros(length(fr), 1);
Ygr = zeros(length(fr), 1);
Ysm = zeros(length(fr), 1);
B = sqrt(w*(w+rs+rr)/((w+rs)*(w+rr)));
for i=1:length(fr)
    waveLen = 340/fr(i);
    gamas = sqrt(2*rs*(w+rr)/(waveLen*L));
    gamar = sqrt(2*rr*(w+rs)/(waveLen*L));
    Ys(i) = gamas.*Mvs;
    Yr(i) = gamar.*Mvr;
    if Ys(i) > Yr(i)
        Ygr(i) = Ys(i);   % greater Y
        Ysm(i) = B*Yr(i);   % smaller Y
    else
        Ygr(i) = Yr(i);   % greater Y
        Ysm(i) = B*Ys(i);   % smaller Y
    end
end
end

 
function d = dist2D(p1, p2)
d = sqrt((p1(1)-p2(1)).^2+(p1(2)-p2(2)).^2);
end