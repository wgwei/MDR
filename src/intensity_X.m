 %load y_cluster stepYgr distributeYgr stepYsm distributeYsm

 if exist('yClusterRott.mat', 'file')~=2
    % load y_cluster for Rotterdam
    %Ygr
    step = 0.1; %the step is fixed as the same as Gent data
    fnsGr = {'srcCateg1_Ygr.txt', 'srcCateg2_Ygr.txt', 'srcCateg3_Ygr.txt', 'srcCateg4_Ygr.txt'};
    fnsSm = {'srcCateg1_Ysm.txt', 'srcCateg2_Ysm.txt', 'srcCateg3_Ysm.txt', 'srcCateg4_Ysm.txt'};

    for g=1:2
        if g==1
            fns = fnsGr;
        else
            fns = fnsSm;
        end
        for f=1:length(fns)
            disp(['processing ' fns{f}]);
            Y = load(fns{f});
            Y = sort(Y);
            distributeY = zeros(length(0:step:Y(end)), length(fns));
            cnt = 1;
            for n=0:step:Y(end)
                distributeY(cnt, f) = length(find(Y<=n & Y>n-step));
                cnt = cnt+1;
            end
        end
        distributeY = sum(distributeY, 2);
        stepY = (0:length(distributeY)-1)*step;
        if g==1
            distributeYgr = distributeY;
            stepYgr = stepY;
        else
            distributeYsm = distributeY;
            stepYsm = stepY;
        end
    end
    save yClusterRott stepYgr distributeYgr stepYsm distributeYsm
 else
     rott = load ('yClusterRott.mat', 'stepYgr', 'distributeYgr', 'stepYsm', 'distributeYsm');
     gent = load ('yClusterGent.mat', 'stepYgr', 'distributeYgr', 'stepYsm', 'distributeYsm');
     
     % check the step
     stepGent = gent.stepYgr(2) - gent.stepYgr(1)
     stepRott = rott.stepYgr(2) - rott.stepYgr(1)
     if stepGent~=stepRott
         disp ('\n\n\t\tThe step is NOT the same!!!')
     else
         dYgr = zeros(max(length(rott.distributeYgr), length(gent.distributeYgr)), 2);
         dYgr(1:length(rott.distributeYgr), 1) = rott.distributeYgr;
         dYgr(1:length(gent.distributeYgr), 2) = gent.distributeYgr;
         dtYgr = sum(dYgr, 2);
         
         dYsm = zeros(max(length(rott.distributeYsm), length(gent.distributeYsm)), 2);
         dYsm(1:length(rott.distributeYsm), 1) = rott.distributeYsm;
         dYsm(1:length(gent.distributeYsm), 2) = gent.distributeYsm;
         dtYsm = sum(dYsm, 2);
         
         % plot out
         intvYgr = (0:length(dtYgr)-1)*stepGent;
         intvYsm = (0:length(dtYsm)-1)*stepGent;
         plot(intvYgr, 100*dtYgr/sum(dtYgr), 'k--', 'linewidth', 1);hold on
         plot(intvYsm, 100*dtYsm/sum(dtYsm), 'k-', 'linewidth', 1);
         legend('Xi+', 'Xi-')
         xlim([0 15]);ylim([0 8])
         xlabel('Input variable');ylabel('Percentage %')
     end          
 end


