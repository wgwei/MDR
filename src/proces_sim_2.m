clear all;
close all;
clc;

addpath('S:\software\intern\script_wgw\special')
addpath('S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\src')
addpath('S:\software\intern\script_wgw\Multiple_diffraction')
dire = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims';


% diffrent comments for different cases.

sfs = { 'BBB-6-12-10-smp64000-wh', ...
        'BBB-6-12-10-smp64000-w',...
        'BTB-10-0-6-smp64000-w',...
        'BB-10-6-smp64000-w'};
sposs = {[5 1], [5 1], [7 1], [7 1]};
vertexs = {[10 24; 16 24; 26 12 ; 33 12; 41 10; 47 10], [10 6; 16 6; 26 12 ; 33 12; 41 10; 47 10],...
    [10 10; 16 10; 20 14 ; 24 6; 30 6], [10 10; 16 10; 24 6; 30 6]};
betas = {[1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi], [1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi], ...
    [1.5*pi 1.5*pi 2*pi 1.5*pi 1.5*pi], [1.5*pi 1.5*pi 1.5*pi 1.5*pi]};
properClaims = {[1 1 1 1 1 1], [1 1 1 1 1 1], [1 1 0 1 1], [1 1 1 1]};
% frs = [31.5 63.5 125 250 500 1000 2000 4000 8000];
frs = [250 2000];
REFDIST = 50;

% plot contour for s=1:length(frs)
for s=1:length(frs)
    fr = frs(s);
    for f=1:length(sfs)
        if f==3 %  'BTB_10_0_6_smp64000_w' -> no adjacet double wedges
            M = 0;
        else
            M = 1;
        end
        sf = sfs{f};
        if exist([dire '\' sf '\' 'lvl_countour_' num2str(fr) 'Hz.mat'], 'file')==2
            load(['lvl_countour_' num2str(fr) 'Hz.mat'], 'calSim', 'surfx', 'surfy')
        else
            spos = sposs{f};
            vertex = vertexs{f};
            beta = betas{f};
            properClaim = properClaims{f};

            cd ([dire '\' sf])    
            height = vertex(end, 2);        
            rp = load ([dire '\' sf '\' ['xverloop2' '.positions.txt']]);    
            rpix = interp1(rp(:,1), rp(:,1), (rp(1,1):0.2:rp(end,1))', 'linear'); % interpolate more receiver position
            surfx = rpix;
            surfy = height:-0.2:0;            

            calSim = zeros(length(surfy), length(surfx));        
            for xv=1:length(surfy)
                rpy = surfy(xv)
                rpi = [rpix, ones(length(rpix), 1)*rpy]; %#ok<AGROW>
                for r= 1:size(rpi,1)
                    rpos = rpi(r, :);
                    srdist = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);    
                    [levelMinusFFatL2, distanceL2] = multi_diffr2_engr(spos, rpos, vertex, beta, properClaim, fr, M, 0);
                    [levelMinusFFatLTheo, distanceLTheo] = multi_diffr2_theo(spos, rpos, vertex, beta, properClaim, fr, M, 0);
                    calSim(xv, r) = levelMinusFFatL2-levelMinusFFatLTheo;
                end
            end
            save(['lvl_countour_' num2str(fr) 'Hz.mat'], 'calSim', 'surfx', 'surfy')               
        end
        figure
        [C,h1] = contour(surfx, surfy, calSim, '-', 'linewidth', 1.5); hold on;
        xlabel('x-coordinate', 'fontsize', 13); ylabel('y-coordinate', 'fontsize', 13)
        set(h1,'ShowText','on','TextStep',get(h1,'LevelStep')*1)
        text_handle = clabel(C, h1);
        set(text_handle, 'fontsize', 13); set(gca,'fontsize', 13)
        title([num2str(fr) 'Hz ' sf])
        dim1 = surfx(end)-surfx(1);
        ylim([surfy(end) surfy(end)+dim1])      
    end
end

