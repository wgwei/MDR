function [theoLevel, simpLevel] = single_wedge_level(nupi, spos, rpos, edgepos, theta, thetas, lamb, varargin)
% calculate the total sound pressure over a knife edge
% based on point source
% nupi is the outside angle of the edge
% thetas is the angle from the edge to the connecting line from source to
%       the diffraction point
% theta is the angle from the edge to the connecting line from receiver to
%       the diffraction point

theta0 = nupi-thetas;
anglePlus = theta+theta0;
angleMinus = theta-theta0;
theoLevel = single_wedge(spos, rpos, edgepos, nupi, anglePlus, angleMinus, lamb);
simpLevel = single_wedge_simplified(spos, rpos, edgepos, nupi, anglePlus, angleMinus, lamb);
end

function vpDirect = p_direct(R, lamb)
    k = 2*pi/lamb;
    vpDirect = exp(1i*k*R)/R;
end

function pdiffLevel = single_wedge(spos, rpos, edgepos, nupi, anglePlus, angleMinus, lamb)   
    r = dist2D(rpos, edgepos);
    r0 = dist2D(spos, edgepos);
    Xplus = X(r, r0, nupi, anglePlus, lamb);
    Xminus = X(r, r0, nupi, angleMinus, lamb);
    pdiffLevel = 10*log10(0.5*((f(Xplus)+f(Xminus))^2+(g(Xplus)+g(Xminus))^2));    
end

function pdiffLevel = single_wedge_simplified(spos, rpos, edgepos, nupi, anglePlus, angleMinus, lamb)   
    r = dist2D(rpos, edgepos);
    r0 = dist2D(spos, edgepos);
    Xplus = X(r, r0, nupi, anglePlus, lamb);
    Xminus = X(r, r0, nupi, angleMinus, lamb);
    pdiffLevel = 10*log10(0.5*((f_simple(Xplus)+f_simple(Xminus))^2));    
end

function vpDiff = p_diff(spos, rpos, edgepos, nupi, anglePlus, angleMinus, lamb)   
    r = dist2D(rpos, edgepos);
    r0 = dist2D(spos, edgepos);
    R = dist2D(spos, rpos);    
    Xplus = X(r, r0, nupi, anglePlus, lamb);
    Xminus = X(r, r0, nupi, angleMinus, lamb);
    ADplus = AD(Xplus);
    ADminus = AD(Xminus);
    L = r+r0;
    k = 2*pi/lamb;
    if abs(ADplus)/abs(ADminus)>5 || abs(ADminus)/abs(ADplus)>5
        yb = ybdline(spos, edgepos, rpos(1));
        N = (L-R)/(0.5*lamb);
        vpDiff = sign(yb-rpos(2))*exp(1i*k*L)/L*exp(1i*pi/4)/sqrt(2)*(f(sqrt(2*N))-1i*g(sqrt(2*N)));
    else
        vpDiff = exp(1i*k*L)/L*exp(1i*pi/4)/sqrt(2)*(ADplus + ADminus);
    end
end

function y = ybdline(spos, edgepos, x)
% calculate the y-coordinate of the boundary line
    y = spos(2) + (edgepos(2)-spos(2))/(edgepos(1)-spos(1))*(x-spos(1));
end

function vAD = AD(X)
    vAD = sign(X)*(f(abs(X))-1i*g(abs(X)));
end

function vx = X(r, r0, nupi, angl, lamb)
    tau = Tau(r,r0, lamb);
    mnv = Mnv(nupi, angl);
    vx = tau*mnv;
end

function mnv = Mnv(nupi, angl)
    nu = pi/nupi;
    mnv = abs((cos(nu*pi)-cos(nu*angl))/(nu*sin(nu*pi)));
end

function tau = Tau(r,r0, lamb)
    tau = sqrt(2*r*r0/(lamb*(r+r0)));
end

function vf = f_simple(x)
    vf = 0.37/(0.37+x);
end

function vf = f(x)
    [c, s] = fcs(x);
    vf  = (0.5-s)*cos(0.5*pi*x^2)-(0.5-c)*sin(0.5*pi*x^2);
end

function vg = g(x)
    [c, s] = fcs(x);
    vg = (0.5-c)*cos(0.5*pi*x^2)+(0.5-s)*sin(0.5*pi*x^2);
end

function d = dist2D(p1, p2)
    d = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
end