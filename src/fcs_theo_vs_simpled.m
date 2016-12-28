% compare theory and simplified diffraction function
clc;close all;clear all;
addpath('S:\software\intern\script_wgw\special')
nv = 2/3;

Ws = 0.01:0.01:10;
dftheo = zeros(1, length(Ws));
dfsimpd = zeros(1, length(Ws));

for n=1:length(Ws)
    w = Ws(n);
    xi = sqrt(0.5*w)*(cos(nv*pi)+1)/nv*sin(nv*pi);
    x2 = sqrt(2/3*w)*(cos(nv*pi)+1)/nv*sin(nv*pi);

    % theory
    [C, S] = fcs(xi);
    fi = (0.5-S)*cos(0.5*pi*xi^2)-(0.5-C)*sin(0.5*pi*xi^2);
    [C2, S2] = fcs(x2);
    f2 = (0.5-S2)*cos(0.5*pi*x2^2)-(0.5-C2)*sin(0.5*pi*x2^2);
    dftheo(n) = fi*f2;

    % simplified
    dfsimpd(n) = 0.37/(0.37+xi)*0.37/(0.37+x2);
end
plot(Ws, 10*log10((dftheo./dfsimpd).^2))


% plot fresnel results
xx=0:0.001:2.5;
Cx = zeros(length(xx), 3);
Sx = zeros(length(xx), 3);
for r=1:length(xx)
    x = xx(r);    
    
    % theory
    [C, S] = fcs(x);
    Cx(r, 1) = C;
    Sx(r, 1) = S;

    % simplified
    Cx(r, 2) = 0.5+0.37/(0.37+x)*sin(0.5*pi*x^2);
    Sx(r, 2) = 0.5-0.37/(0.37+x)*cos(0.5*pi*x^2);
    
    % simplified 2
    Cx(r, 3) = 0.5+1/(pi*x)*sin(0.5*pi*x^2);
    Sx(r, 3) = 0.5-1/(pi*x)*cos(0.5*pi*x^2);
end

figi = figure;
set(figi,'units','inches', 'Position', [2, 2, 11, 4])

subplot(1, 2, 1)
sybs = {'-', '--', '-.'};
for p=1:3
    plot(xx,Cx(:, p), sybs{p}, 'linewidth', 1); hold on
end
legend('Theory', '0.5+0.37/(x+0.37)*sin(0.5x^2)',  '0.5+1/(pi*x)*sin(0.5x^2)', 'location', 'best')
xlabel('x')
ylabel('Fresnel integral C(x)')
ylim([0 1])

subplot(1, 2, 2)
for p=1:3
    plot(xx, Sx(:,p), sybs{p}, 'linewidth', 1); hold on
end
legend('Theory', '0.5-0.37/(x+0.37)*cos(0.5x^2)',  '0.5-1/(pi*x)*cos(0.5x^2)', 'location', 'best')
xlabel('x')
ylabel('Fresnel integral S(x)')
ylim([-1,1])


% plot the diffraction function
