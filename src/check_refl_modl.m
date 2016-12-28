% check the calculation of the reflection method
clear all; 
close all;
clc;

addpath('S:\software\intern\script_wgw\special');
addpath('S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\script')
addpath('S:\software\intern\script_wgw\double_diffraction')
ratio = (1:2:49)./(3:2:51);
fr = [31.5 63 125 250 500 1000];
alfa = 0.97;
spos = [10, 1];
rs = sqrt(10^2+11^2);
thetas = atan(11/10);
Ws = 20;

% for Hs==Hi case
Hs = 12;
h1 = 11;
Hi = 12;
W23 = 10;
edgePos = [30 12];
dire = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BB_12_12smp64000';
xvs = 2:12;
for am=1:length(fr)
    yy = xvs-1+0.81;
    rcvFile = ['xverloop' num2str(2) '.positions.txt'];
    rposs = load([dire '\' rcvFile]);
    xx = rposs(:, 1);
    z2 = zeros(length(yy), length(xx));
    for v=1:length(xvs)
        xv  = xvs(v);
        rcvFile = ['xverloop' num2str(xv) '.positions.txt'];
        vFDTD = load ([dire '\' 'xverloop'  num2str(xv) '_L2P_oct.mat'], 'Lrf_oct');
        rposs = load([dire '\' rcvFile]);
        for r=1:size(rposs, 1)
            rpos = rposs(r, :);
            rr = sqrt((edgePos(1)-rpos(1))^2 + (edgePos(2)-rpos(2))^2);
            srdist = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);
            phir = atan(abs((rpos(1)-edgePos(1))/(edgePos(2)-rpos(2))));        
            lamb = 340/fr(am);            
            
            C1s = (0.37./(sqrt(3)*0.5*sqrt(2*rr./lamb).*cos(phir)+0.37)).^2;
            C3s = 3.31*h1*sqrt(W23/lamb)+0.5*Ws+rr+W23;
            psq = C1s.*alfa^2/Ws^2 *lerch(alfa^2, 2, (C3s+Ws)/Ws);
            
            intnRatio = intn_over_ff_double(rs, rr, W23, 1.5*pi, 1.5*pi, thetas, phir, lamb);
            z2(v,r) = 10*log10((psq + intnRatio/(rs+rr+W23)^2 )*srdist^2)-...
                vFDTD.Lrf_oct(am, r);
        end        
    end
    figure
    [C, hd1] = contour(xx, yy, z2);
    set(hd1,'ShowText','on','TextStep',get(hd1,'LevelStep')*1)
%     xlabel('Horizontal Location m', 'fontsize', 13); ylabel('Vertical Location m', 'fontsize', 13)
    xlabel('Horizontal Location m'); ylabel('Vertical Location m')
    title([num2str(fr(am)) ' Hz'])
%     clabel(C, hd1, 'fontsize', 12)
%     set(gca,  'fontsize', 13)
    saveas(gcf, ['fig\HH_' num2str(fr(am)) 'Hz.fig'])
end


% for Hs<Hi case
Hs = 8;
h1 = 11;
Hi = 12;
W23 = 10;
edgePos = [30 12];
epos = [20 12];
dire = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BB_8_12_smp64000';
xvs = 2:12;
for am=1:length(fr)
    yy = xvs-1+0.81;
    rcvFile = ['xverloop' num2str(2) '.positions.txt'];
    rposs = load([dire '\' rcvFile]);
    xx = rposs(:, 1);
    z2 = zeros(length(yy), length(xx));
    for v=1:length(xvs)
        xv  = xvs(v);
        rcvFile = ['xverloop' num2str(xv) '.positions.txt'];
        vFDTD = load ([dire '\' 'xverloop'  num2str(xv) '_L2P_oct.mat'], 'Lrf_oct');
        rposs = load([dire '\' rcvFile]);
        for r=1:size(rposs, 1)
            rpos = rposs(r, :);
            rr = sqrt((edgePos(1)-rpos(1))^2 + (edgePos(2)-rpos(2))^2);
            srdist = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);
            phir = atan(abs((rpos(1)-edgePos(1))/(edgePos(2)-rpos(2))));        
            lamb = 340/fr(am);
            W12 = sqrt(Ws^2+(Hi-Hs).^2);
            psi = acos((Hi-Hs)./W12);
            C1s = (0.37./(sqrt(3)*0.5*sqrt(2*rr./lamb).*cos(phir)+0.37)).^2;
            C1sp = (0.37./(2/3.*sqrt(2*W12./lamb).*cos(psi)+0.37)).^2 .* ...
                (0.37./(sqrt(3)*0.5*sqrt(2*W23./lamb).*cos(psi)+0.37)).^2 .* ...
                (0.37./(sqrt(3)*0.5*sqrt(2*rr./lamb).*cos(phir)+0.37)).^2;
            C3s = 0.5*Ws+W23+rr+sqrt(2*W23/lamb)*sqrt(3)/0.74*h1;
            C3sp = 0.5*Ws+W23+rr+W12;
            idx = find((Hs-spos(2))/(Hi-spos(2))<ratio);
            if length(idx)>=1
                N = idx(1);
            else
                N = length(ratio);
            end  
            intnRatio = intn_over_ff_double(rs, rr, W23, 1.5*pi, 1.5*pi, thetas, phir, lamb);
            z2(v,r) = 10*log10((C1sp*alfa^2/Ws^2*lerch(alfa^2, 2, 1+C3sp/Ws) + ...
                intnRatio/(rs+rr+W23)^2 + ...
                C1s*alfa^2/Ws^2*lerch(alfa^2, 2, 1+C3s/Ws) - ...
                C1s*alfa^(2*(N+1))/Ws^2*lerch(alfa^2, 2, N+1+C3s/Ws))*srdist^2)-...
                vFDTD.Lrf_oct(am, r);
        end        
    end
    figure
    [C, hd1] = contour(xx, yy, z2);
    set(hd1,'ShowText','on','TextStep',get(hd1,'LevelStep')*1)
%     xlabel('Horizontal Location m', 'fontsize', 13); ylabel('Vertical Location m', 'fontsize', 13)
    xlabel('Horizontal Location m'); ylabel('Vertical Location m')
    title([num2str(fr(am)) ' Hz'])
%     clabel(C, hd1, 'fontsize', 12)
%     set(gca,  'fontsize', 13)
    saveas(gcf, ['fig\LH_' num2str(fr(am)) 'Hz.fig'])
end


% for Hs>Hi case
alfa = 0.97;
Ws2 = 20;
Hs2 = 12;
h1 = 7;
Hi2 = 8;
Wi = 10;
rs2 = sqrt(10^2+7^2);
thetas2 = atan(7/10);
z2 = zeros(length(fr), 1);
edgePos2 = [30 8];
epos = [20 8];
dire2 = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BB_12_8_smp64000';
xvs = 2:8;
omiga = atan(abs(Hi2-Hs2)/(Ws2+Wi));
for am=1:length(fr)
    yy = xvs-1+0.81;
    rcvFile = ['xverloop' num2str(2) '.positions.txt'];
    rposs = load([dire2 '\' rcvFile]);
    xx = rposs(:, 1);
    z2 = zeros(length(yy), length(xx));
    for v=1:length(xvs)
        xv  = xvs(v);
        rcvFile = ['xverloop' num2str(xv) '.positions.txt'];
        vFDTD = load ([dire2 '\' 'xverloop'  num2str(xv) '_L2P_oct.mat'], 'Lrf_oct');
        rposs = load([dire2 '\' rcvFile]);
        for r=1:size(rposs, 1)
            rpos = rposs(r, :);
            rr = sqrt((edgePos2(1)-rpos(1))^2 + (edgePos2(2)-rpos(2))^2);
            srdist = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);
            phir = atan(abs((rpos(1)-edgePos2(1))/(edgePos2(2)-rpos(2))));        
            lamb = 340/fr(am);
            W12_2 = sqrt((Ws2+Wi)^2+(Hi2-Hs2).^2);
            psi = acos((Hi2-Hs2)./W12_2);
            C1s = (0.37./(sqrt(3)*0.5*sqrt(2*rr./lamb).*cos(phir)+0.37)).^2;
            C1 = (0.37./(sqrt(2*3*rr/lamb)*(cos(2/3*phir)-0.5)+0.37)).^2;
            C2 = abs(Ws2-1.80*Ws2*sqrt(2*W12_2/lamb)*sin(omiga));
            C3 = 1.80*sqrt(2*W12_2/lamb)*(Hs2-spos(2))*cos(omiga)-1.80*sqrt(2*W12_2/lamb)*0.5*Ws2*sin(omiga)+0.5*Ws2+W12_2+rr;
            C3s = 0.5*Ws2+Wi+rr+sqrt(2*Wi/lamb)*sqrt(3)/0.74*h1;
            C3sp = 0.5*Ws2+Wi+rr+W12_2;
            idx = find((Hi2-spos(2))/(Hs2-spos(2))<ratio);
            
            alfa5 = atan(abs(Hi2-Hs2)/(Ws2+Wi));
            alfa6 = 1.5*pi - phir;
            alfa3 = pi + psi;
            alfa4 = pi + alfa5;
            CNT1 = (0.37/(sqrt(W12_2*2*3/lamb)*abs(0.5+cos(alfa4-alfa3)) + 0.37))^2 ...
            * (0.37/(sqrt(rr*2*3/lamb)*abs(0.5+cos(alfa6-alfa5)) + 0.37))^2;
            CNT2 = 1.8*sqrt(2*W12/lamb)*sin(psi)*Ws2+Ws2;
            CNT3 = 1.8*sqrt(2*W12/lamb)*(h1*cos(psi)+0.5*Ws2*sin(psi));
            if length(idx)>=1
                N = idx(1);
            else
                N = length(ratio);
            end
            intnRatio = intn_over_ff_double(rs2, rr, W23, 1.5*pi, 1.5*pi, thetas2, phir, lamb);
            z2(v,r) = 10*log10((C1*alfa^2/C2^2*lerch(alfa^2, 2, 1+C3/C2) - ...   
                C1*alfa^(2*(N+1))/C2^2*lerch(alfa^2, 2, N+1+C3/C2) + ...
                C1s*alfa^2/Ws2^2*lerch(alfa^2, 2, 1+C3s/Ws2) + ...
                intnRatio/(rs2+rr+Wi)^2 + ...
                CNT1*alfa^2/CNT2^2*lerch(alfa^2, 2, 1+CNT3/CNT2))*srdist^2)-...
                vFDTD.Lrf_oct(am, r);
        end        
    end
    figure
    [C, hd1] = contour(xx, yy, z2);
    set(hd1,'ShowText','on','TextStep',get(hd1,'LevelStep')*1)
%     xlabel('Horizontal Location m', 'fontsize', 13); ylabel('Vertical Location m', 'fontsize', 13)
    xlabel('Horizontal Location m'); ylabel('Vertical Location m')
%     title([num2str(fr(am)) ' Hz'])
%     clabel(C, hd1, 'fontsize', 12)
%     set(gca,  'fontsize', 13)
    saveas(gcf, ['fig\HL_' num2str(fr(am)) 'Hz.fig'])
end

% multiple reflections and diffraction
Hs = 12;
h1 = 11;
Hi = 12;
W12 = 10.;
W23 = sqrt(256+36);
edgePos = [46 6];
epos = [30 12];
dire = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BBB-12-12-6-smp64000';
xvs = 2:10;
for am=1:length(fr)
    lamb = 340/fr(am);
    yy = [1.81 2.81 3.81 4.81 5.81 6.81 7.81 8.81 9.81 10.81 11.91];
    rcvFile = ['xverloop' num2str(2) '.positions.txt'];
    rposs = load([dire '\' rcvFile]);
    xx = rposs(:, 1);
    z2 = zeros(length(yy), length(xx));
    for v=1:length(xvs)
        xv  = xvs(v);
        rcvFile = ['xverloop' num2str(xv) '.positions.txt'];
        vFDTD = load ([dire '\' 'xverloop'  num2str(xv) '_L2P_oct.mat'], 'Lrf_oct');
        rposs = load([dire '\' rcvFile]);
        for r=1:size(rposs, 1)
            rpos = rposs(r, :) ;     
            srdist = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);
            [ix,iy] = perpendicularPnt(epos, rpos, edgePos); % perpendicular intersection
            d1 = sqrt((epos(1)-ix)^2+(epos(2)-iy)^2);
            d2 = sqrt((rpos(1)-ix)^2+(rpos(2)-iy)^2);
            d3 = sqrt((edgePos(1)-ix)^2+(edgePos(2)-iy)^2);
            rfn = sqrt(lamb*d1*d2/(d1+d2)); % radius of fresnel zone
            if iy>=edgePos(2) % double diffraction
                rr = sqrt((epos(1)-rpos(1))^2 + (epos(2)-rpos(2))^2);
                phir = atan(abs((rpos(1)-epos(1))/(epos(2)-rpos(2)))); 
                C1s = (0.37./(sqrt(3)*0.5*sqrt(2*rr./lamb).*cos(phir)+0.37)).^2;
                C3s = 0.5*Ws+W12+rr+3.31*sqrt(W12/lamb)*h1;         
                if d3>rfn
                    fresnelCor = 0;
                else
                    [sperc, unsperc] = fresnelperc(lamb, epos, rpos, edgePos); % only for this special case 
                    fresnelCor = 10*log10(unsperc);
                end
                intnRatio = intn_over_ff_double(rs, rr, W12, 1.5*pi, 1.5*pi, thetas, phir, lamb);
                z2(v,r) = 10*log10((C1s*alfa^2/Ws^2*lerch(alfa^2, 2, 1+C3s/Ws) + ...
                    intnRatio/(rs+rr+W12)^2)*srdist^2) - vFDTD.Lrf_oct(am, r)+fresnelCor;
            else % 3-order diffraction
                rr = sqrt((edgePos(1)-rpos(1))^2 + (edgePos(2)-rpos(2))^2);
                phir = atan(abs((rpos(1)-edgePos(1))/(edgePos(2)-rpos(2))));  
                theta2 = atan(16/6); % angle between two building connecting line # it's only for this test case
                theta3 = 1.5*pi-atan(6/16); % ouside angle between two building connecting line # it's only for this test case
                nu = 2/3;
                Mvplus = (cos(nu*pi)-cos(nu*(1.5*pi+theta2)))/(nu*sin(nu*pi));
                Mvminus = (cos(nu*pi)-cos(nu*(1.5*pi-theta2)))/(nu*sin(nu*pi));
                MvplusR = (cos(nu*pi)-cos(nu*(phir+theta3)))/(nu*sin(nu*pi));
                MvminusR = (cos(nu*pi)-cos(nu*(phir-theta3)))/(nu*sin(nu*pi));
                C1s = 0.25*(0.37./(sqrt(2*W23./lamb).*Mvplus+0.37)+0.37./(sqrt(2*W23./lamb).*Mvminus+0.37)).^2 * ...
                    (0.37./(sqrt(2*rr./lamb).*MvplusR+0.37)+0.37./(sqrt(2*rr./lamb).*MvminusR+0.37)).^2;
                C3s = 0.5*Ws+W12+W23+rr+3.31*sqrt(W12/lamb)*h1;
                if d3>rfn
                    fresnelCor = 0;
                else
                    [sperc, unsperc] = fresnelperc(lamb, epos, rpos, edgePos); % only for this special case 
                    fresnelCor = 10*log10(sperc);
                end
                intnRatio = multiple_dif_fun_v3(rs, rr, [W12 W23], lamb, ...
                    [1.5*pi-atan(10/11) 1.5*pi 1.5*pi-atan(6/16)], [0 atan(16/6) phir], [1.5*pi 1.5*pi 1.5*pi]);
                z2(v,r) = 10*log10((C1s*alfa^2/Ws^2*lerch(alfa^2, 2, 1+C3s/Ws) + ...
                    0.25*(abs(intnRatio).^2)./(rs+rr+W12+W23)^2)*srdist^2) - vFDTD.Lrf_oct(am, r)+fresnelCor; % 0.25 = (1/2)^C where C=2
            end            
        end        
    end
    figure    
    [C, hd1] = contour(xx, yy, z2);
    set(hd1,'ShowText','on','TextStep',get(hd1,'LevelStep')*1)
%     xlabel('Horizontal Location m', 'fontsize', 13); ylabel('Vertical Location m', 'fontsize', 13)
    xlabel('Horizontal Location m'); ylabel('Vertical Location m')
%     title([num2str(fr(am)) ' Hz'])
%     clabel(C, hd1, 'fontsize', 12)
%     set(gca,  'fontsize', 13)
    saveas(gcf, ['fig\COMB_' num2str(fr(am)) 'Hz.fig'])
end