% compare theory and simplified diffraction function
clc;close all;clear all;

% single diffraction
rs = 10;
rr = 10;
spos = [-10*sin(pi/3), -10*sin(pi/6)];
rpos = [10 0]; % in this comparison rpos is only used to calculate r and r0. so it's value can be fixed
edgepos = [0 0];
nupi = 11/6*pi;
thetas = pi/6;
Theta = 0:0.1:2*pi/3;
theoLevel = zeros(1, length(Theta));
simpLevel = zeros(1, length(Theta));
figi = figure;
lambs = [0.1 1 10];
for n=1:3
    lamb = lambs(n);
    for t=1:length(Theta)
        theta = Theta(t);
        [tv, sv] = single_wedge_level(nupi, spos, rpos, edgepos, theta, thetas, lamb);
        theoLevel(t) = tv;
        simpLevel(t) = sv;
    end
    set(figi,'units','inches', 'Position', [2, 2, 11, 4])
    subplot(1, 2, 1)
    plot(Theta, theoLevel, '-', 'linewidth', 1);hold on
    plot(Theta, simpLevel, '--', 'linewidth', 1)
    xlabel('Theta [radian]')
    ylabel('Relative sound pressure level [dB]')
end
legend('Pierce', 'Simplified', 'location', 'best')

% double diffraction
spos = [-10*sqrt(2)/2, -10*sqrt(2)/2];
N1p = [0 0];
N2p = [10 0];
betas = 1.5*pi;
betar = 1.5*pi;

thetar = 0:0.01:pi/2;
lvlSimp = zeros(length(thetar), 1);
lvlTheo = zeros(length(thetar), 1);

frs = 340./[10 1 0.1];

for c=1:3
    fr = frs(c);
    for n=1:length(thetar)
        ag = thetar(n);
        rpos = [10+10*sin(ag) -10*cos(ag)];
        [p2lvlSimpd, p2lvlTheo, L] = double_diffr(spos, rpos, N1p, N2p, fr, betas, betar);
        lvlSimp(n) = p2lvlSimpd;
        lvlTheo(n) = p2lvlTheo;
    end
    subplot(1, 2, 2)
    plot(thetar, lvlSimp, '--', 'linewidth', 1); hold on
    plot(thetar, lvlTheo, '-', 'linewidth', 1)    
    xlabel('Theta_r [radian]')
    ylabel('Relative sound pressure level [dB]')
    v = axis;
end
legend('Pierce', 'Simplified', 'location', 'best')
% set(handle,'Position',[v(1)*1.2 v(4)*.8 0]); 
