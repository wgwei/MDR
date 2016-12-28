classdef error_imag_src 
% calculate the contribution of image sources
% plot out the results: error vs number of image sources
    properties
        nu1 = 2/3;
        nu2 = 2/3;
        beta1 = 1.5*pi;
        beta2 = 1.5*pi;
        phir = pi/4;
        rr = 7;
        w = 22;
        lam = 0.68;
        alfa = 0.9744;
    end
    methods 
        function obj = error_imag_src()
            addpath('S:\software\intern\script_wgw\special')
%             cd('S:\projecten\QSIDE\Action 2\background mapping\submitted_publications\Wei_engineeringModel\scripts\compare_image_source');
            call_plot(obj)
        end
        
        function [Ldp_n, Ldp] = L_imag50(obj, Nr, ws, h1, alfa)
            % sum up contributioin of Nr image sources
            % Ldp_n: level ref L=rs+w+rr between source, bar, receiver;
            % Ldp: is the final value of Nr image sources.
            psq = zeros(Nr, 1);
            L = sqrt((0.5*ws)^2+h1^2)+obj.w+obj.rr;
            for n =1:Nr
                rs = sqrt((n*ws+0.5*ws)^2+h1^2);
                phis = atan((n*ws+0.5*ws)/h1);
                gammas = sqrt(2*rs*(obj.w+obj.rr)/(obj.lam*(rs+obj.rr+obj.w)));
                Mvs = (cos(obj.nu1*pi)-cos(obj.nu1*(obj.beta1-phis)))/(obj.nu1*sin(obj.nu1*pi));
                gammar = sqrt(2*obj.rr*(obj.w+rs)/(obj.lam*(rs+obj.rr+obj.w)));
                Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-obj.phir)))/(obj.nu2*sin(obj.nu2*pi));
                B = sqrt(obj.w*(obj.w+obj.rr+rs)/((obj.w+rs)*(obj.w+obj.rr)));
                Ys = gammas*Mvs;
                Yr = gammar*Mvr;
                if Ys<Yr
                    Ysm = B*Ys;
                    Ygr = Yr;
                else
                    Ysm = B*Yr;
                    Ygr = Ys;
                end
                Ln = rs+obj.w+obj.rr;
                psq(n) = alfa^(n*2) * (0.37/(Ysm+0.37))^2 * (0.37/(Ygr+0.37))^2 / Ln^2;     
            end
            Ldp_n = 10*log10(cumsum(psq)* L^2);
            Ldp = Ldp_n(end);
        end        
            
        function [Ldp_n, Ldp] = L_imag50byInput(obj, Nr, ws, w, rr, h1, phir, lam, alfa)
            % sum up contributioin of Nr image sources
            % Ldp_n: level ref L=rs+w+rr between source, bar, receiver;
            % Ldp: is the final value of Nr image sources.
            psq = zeros(Nr, 1);
            L = sqrt((0.5*ws)^2+h1^2)+w+rr;
            for n =1:Nr
                rs = sqrt((n*ws+0.5*ws)^2+h1^2);
                phis = atan((n*ws+0.5*ws)/h1);
                gammas = sqrt(2*rs*(w+rr)/(lam*(rs+rr+w)));
                Mvs = (cos(obj.nu1*pi)-cos(obj.nu1*(obj.beta1-phis)))/(obj.nu1*sin(obj.nu1*pi));
                gammar = sqrt(2*rr*(w+rs)/(lam*(rs+rr+w)));
                Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
                B = sqrt(w*(w+rr+rs)/((w+rs)*(w+rr)));
                Ys = gammas*Mvs;
                Yr = gammar*Mvr;
                if Ys<Yr
                    Ysm = B*Ys;
                    Ygr = Yr;
                else
                    Ysm = B*Yr;
                    Ygr = Ys;
                end
                Ln = rs+w+rr;
                psq(n) = alfa^(n*2) * (0.37/(Ysm+0.37))^2 * (0.37/(Ygr+0.37))^2 / Ln^2;     
            end
            Ldp_n = 10*log10(cumsum(psq)* L^2);
            Ldp = Ldp_n(end);
        end        
        
        function L1to = L_imag_approx(obj, ws, h1, alfa)
            L = sqrt((0.5*ws)^2+h1^2)+obj.w+obj.rr;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-obj.phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*obj.rr/obj.lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C2s = 0.5*ws+obj.w+obj.rr;
            L1to = 10*log10(C1s*alfa^2/ws^2*lerch(alfa^2, 2, (C2s+ws)/ws)*L^2);
        end
        
        function Ldp = L_approx_refl(obj, ws, h1, alfa)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = sqrt((0.5*ws)^2+h1^2)+obj.w+obj.rr;
            Mvr = sqrt(3)*(cos(2/3*obj.phir)-0.5);
            gammar = sqrt(2*obj.rr/obj.lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(obj.w/obj.lam)+1.5*ws+obj.rr+obj.w;
            psq = C1s*alfa^2/ws^2*lerch(alfa^2, 2, (C3s+ws)/ws);
            Ldp = 10*log10(psq*L^2);
        end
        
        function Ldp00 = Lp00(obj,  ws, h1, alfa)
            rs = sqrt((0.5*ws)^2+h1^2);
            L = rs+obj.w+obj.rr;            
            phis = atan(0.5*ws/h1);
            gammas = sqrt(2*rs*(obj.w+obj.rr)/(obj.lam*L));
            Mvs = (cos(obj.nu1*pi)-cos(obj.nu1*(obj.beta1-phis)))/(obj.nu1*sin(obj.nu1*pi));
            gammar = sqrt(2*obj.rr*(obj.w+rs)/(obj.lam*L));
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-obj.phir)))/(obj.nu2*sin(obj.nu2*pi));
            B = sqrt(obj.w*(obj.w+obj.rr+rs)/((obj.w+rs)*(obj.w+obj.rr)));
            Ys = gammas*Mvs;
            Yr = gammar*Mvr;
            if Ys<Yr
                Ysm = B*Ys;
                Ygr = Yr;
            else
                Ysm = B*Yr;
                Ygr = Ys;
            end
            Ldp00 = 10*log10(alfa^2 * (0.37/(Ysm+0.37))^2 * (0.37/(Ygr+0.37))^2);     
        end
        
        function [LbarMff, Ln] = Lb_ff(obj, ws, w, rr, h1, phir, lam, alfa)
            rs = sqrt((0.5*ws)^2+h1^2);
            phis = atan(0.5*ws/h1);
            gammas = sqrt(2*rs*(w+rr)/(lam*(rs+rr+w)));
            Mvs = (cos(obj.nu1*pi)-cos(obj.nu1*(obj.beta1-phis)))/(obj.nu1*sin(obj.nu1*pi));
            gammar = sqrt(2*rr*(w+rs)/(lam*(rs+rr+w)));
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            B = sqrt(w*(w+rr+rs)/((w+rs)*(w+rr)));
            Ys = gammas*Mvs;
            Yr = gammar*Mvr;
            if Ys<Yr
                Ysm = B*Ys;
                Ygr = Yr;
            else
                Ysm = B*Yr;
                Ygr = Ys;
            end
            Ln = rs+w+rr;
            psq = alfa^2 * (0.37/(Ysm+0.37))^2 * (0.37/(Ygr+0.37))^2;    
            LbarMff = 10*log10(psq);
        end  
        
        function Ldp = L_approx_reflbyInput(obj, alfa, ws, wi, rs, rr, h1, phir, lam)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = rs+rr+wi;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*rr/lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(wi/lam)+1.5*ws+rr+wi;
            psq = C1s*(alfa*sqrt(2)/2)^2/ws^2*lerch((alfa*sqrt(2)/2)^2, 2, (C3s+ws)/ws);
            Ldp = 10*log10(psq*L^2);
        end
        function Ldp = L_approx_reflbyInput1(obj, alfa, ws, wi, rs, rr, h1, phir, lam)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = rs+rr+wi;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*rr/lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(wi/lam)+1.5*ws+rr+wi;
            psq = C1s*alfa^2/ws^2*1.59*alfa^2/((C3s+ws)/ws)^2;
            Ldp = 10*log10(psq*L^2);
        end
        function Ldp = L_approx_reflbyInput2(obj, alfa, ws, wi, rs, rr, h1, phir, lam)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = rs+rr+wi;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*rr/lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(wi/lam)+1.5*ws+rr+wi;
            psq = C1s*alfa^2/ws^2*2.15*alfa^2/(((C3s+ws)/ws)^2+0.13);
            Ldp = 10*log10(psq*L^2);
        end
         function Ldp = L_approx_reflbyInput3(obj, alfa, ws, wi, rs, rr, h1, phir, lam)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = rs+rr+wi;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*rr/lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(wi/lam)+1.5*ws+rr+wi;
            psq = C1s*alfa^2/ws^2*2.72*alfa^2/((C3s+ws)/ws+0.19)^2;
            Ldp = 10*log10(psq*L^2);
         end
         function Ldp = L_approx_reflbyInput4(obj, alfa, ws, wi, rs, rr, h1, phir, lam)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = rs+rr+wi;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*rr/lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(wi/lam)+1.5*ws+rr+wi;
            psq = C1s*alfa^2/ws^2*1.68*alfa^2/(((C3s+ws)/ws)^1.48-0.073);
            Ldp = 10*log10(psq*L^2);
         end
         function Ldp = L_approx_reflbyInput5(obj, alfa, ws, wi, rs, rr, h1, phir, lam)
            % calculate the level reference level in free field with 
            % distance between S->1->2->R. 
            L = rs+rr+wi;
            Mvr = (cos(obj.nu2*pi)-cos(obj.nu2.*(obj.beta2-phir)))/(obj.nu2*sin(obj.nu2*pi));
            gammar = sqrt(2*rr/lam);
            C1s = (0.37/(Mvr*gammar+0.37))^2;
            C3s = 3.31*h1*sqrt(wi/lam)+1.5*ws+rr+wi;
            psq = C1s*alfa^2/ws^2*1.55*alfa^2/((C3s+ws)/ws-0.104)^1.44;
            Ldp = 10*log10(psq*L^2);
         end
        function plot_approx_50imgsrc(obj)
            h1 = 5:20;  % 5
            ws = 6:40;  % 3
            nh1 = length(h1);
            nws = length(ws);
            Lapprox = zeros(nws, nh1);
            Lapprox2 = zeros(nws, nh1);
            Limgsrc50 = zeros(nws, nh1);
            Ldiff = zeros(nws, nh1);       
            Ldiff2 = zeros(nws, nh1);
            for m = 1:nh1
                h1v = h1(m);
                for n = 1:nws
                    wsv = ws(n);
                    Lapprox(n, m) = obj.L_imag_approx(wsv, h1v, obj.alfa);
                    Lapprox2(n, m) = obj.L_approx_refl(wsv, h1v, obj.alfa);
                    [aaa, Limgsrc50(n, m)] = obj.L_imag50(50, wsv, h1v, obj.alfa);
                    Ldiff(n, m) = Lapprox(n, m) - Limgsrc50(n, m);
                    Ldiff2(n, m) = Lapprox2(n, m) - Limgsrc50(n, m);
                end
            end
            figure
            surf(h1, ws, Ldiff)
            xlabel('h1'); ylabel('W_s'); zlabel('Lapprox-Limgsrc50')
            figure
            surf(h1, ws, Ldiff2)
            xlabel('h1'); ylabel('W_s'); zlabel('Lapprox2-Limgsrc50')
        end
        
        function plot_50imgsrc(obj)
            h1 = 10;
            ws = 22;      
            alfs = [0.8 0.9 0.98];
            syms = {'k-', 'b--', 'r-.'};
            subplot(1,2,1)
            for m = 1:3
                alf = alfs(m);
                [Ldp_n, Ldp] = obj.L_imag50(50, ws, h1, alf);                %#ok<NASGU>
                plot(1:50, Ldp_n, syms{m}); hold on;
            end
            grid on; legend('refl coef=0.8', 'refl coef=0.9', 'refl coef=0.98');
            xlabel('Number of image sources');
            ylabel('Level minus free field level at L dB');
        end
        
        function plot_1stimag_2to(obj)
            h1 = 5:20;  % 5
            ws = 6:40;  % 3
            nh1 = length(h1);
            nws = length(ws);
            L2to = zeros(nws, nh1);
            L1 = zeros(nws, nh1);

            for m = 1:nh1
                h1v = h1(m);
                for n = 1:nws
                    wsv = ws(n);
                    [Ldp_n,Ldp] = obj.L_imag50(50, wsv, h1v, obj.alfa);
                    L2to(n, m) = 10*log10(10^(0.1*Ldp) - 10^(0.1*Ldp_n(1)));
                    L1(n, m) = Ldp_n(1);                      
                end
            end
            subplot(1, 2, 2)
            surf(h1, ws, L2to-L1)
            shading interp
            xlabel('h [m]'); ylabel('W_s [m]'); zlabel('Level difference dB')
        end
        
        function plot_cmp_FDTD(obj)
            dire = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BB_12_12smp64000';            
            spos = [10.00 1.00];
            bding = [20 12; 30 12];            
            h1 = bding(1, 2)-spos(2);
            rs = sqrt(h1^2+(bding(1, 1)-spos(1))^2);
            wi = 10;
            ws = 20;            
            frs = [63. 125. 250. 500. 1000.];
            rposen = load ([dire '\' 'xverloop2.positions.txt']);
            [r, aaa] = size(rposen);
            nws = length(frs);
            L1toAppr = zeros(nws, r);
            L1toAppr2 = zeros(nws, r);
            L50imgs = zeros(nws, r);
            for m = 1:r
                rpos = rposen(m, :);
                dst = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);
                rri = sqrt((rpos(1,1)-bding(2,1))^2+(rpos(1,2)-bding(2,2))^2);                
                for n = 1:nws                    
                    fr = frs(n);
                    lami = 340/fr;                    
                    phiri = atan(abs(rpos(1,1)-bding(2,1))/abs(rpos(1,2)-bding(2,2)));
                    [LbarMff, aaa] = Lb_ff(obj, ws, wi, rri, h1, phiri, lami, obj.alfa);
%                     k = 2*pi/lami; 
%                     Lcdst = 10*log10(abs(besselh(0,1,k*dst))/abs(besselh(0,1,k*(rri+wi+rs))));
                    Lcdst2 = 20*log10(dst/(rri+wi+rs));
                    [aaa, Ldp] = obj.L_imag50byInput(50, ws, wi, rri, h1, phiri, lami, obj.alfa);                         
                    Lappr = L_imag_approxbyInput(obj,ws, wi, rri, h1, phiri, lami, obj.alfa);
                    Lappr2 = L_approx_reflbyInput(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);                   
                    
                    % line source should be Lcdst. But the approximation of the reflection is based on point source. 
                    % Then the correction becomes the point source correction, which should be added on the source. 
                    % Here just put the correction up, there is no physical
                    % meaning
                    L1toAppr(n, m) = 10*log10(10^(0.1*(Lappr))+10^(0.1*(LbarMff)))+Lcdst2; 
                    L1toAppr2(n, m) = 10*log10(10^(0.1*(Lappr2))+10^(0.1*(LbarMff)))+Lcdst2;
                    L50imgs(n, m) = 10*log10(10^(0.1*(Ldp))+10^(0.1*(LbarMff)))+Lcdst2; 
                end
            end
            % the following data is based on point source by timing 1/sqrt(ct)
            % to the line source.
            load('S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BB_12_12smp64000\levelminusFFL2P.mat', 'Lrf_oct');
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-L1toAppr) %#ok<NODEF> % start from 63Hz
            xlabel('recvr position'); ylabel('frequency'); zlabel('FDTD-L1toAppr')
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-L1toAppr2)
            xlabel('recvr position'); ylabel('frequency'); zlabel('FDTD-L1toAppr2')
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-L50imgs)
            xlabel('recvr position'); ylabel('frequency'); zlabel('FDTD-L50imgs')
        end       
        
        function plot_cmp_FDTD_lerch(obj)
            dire = 'S:\literatuur\references_wgw\PhDstuff\Publications\Multiple_diffraction\sims\BB_12_12smp64000';            
            spos = [10.00 1.00];
            bding = [20 12; 30 12];            
            h1 = bding(1, 2)-spos(2);
            rs = sqrt(h1^2+(bding(1, 1)-spos(1))^2);
            wi = 10;
            ws = 20;            
            frs = [63. 125. 250. 500. 1000.];
            rposen = load ([dire '\' 'xverloop2.positions.txt']);
            [r, aaa] = size(rposen);
            nws = length(frs);
            approLerch = zeros(nws, r);
            appro1 = zeros(nws, r);            
            appro2 = zeros(nws, r);
            appro3 = zeros(nws, r);
            appro4 = zeros(nws, r);
            appro5 = zeros(nws, r);
            for m = 1:r
                rpos = rposen(m, :);
                dst = sqrt((spos(1)-rpos(1))^2+(spos(2)-rpos(2))^2);
                rri = sqrt((rpos(1,1)-bding(2,1))^2+(rpos(1,2)-bding(2,2))^2);                
                for n = 1:nws                    
                    fr = frs(n);
                    lami = 340/fr;                    
                    phiri = atan(abs(rpos(1,1)-bding(2,1))/abs(rpos(1,2)-bding(2,2)));
                    [LbarMff, aaa] = Lb_ff(obj, ws, wi, rri, h1, phiri, lami, obj.alfa);
%                     k = 2*pi/lami; 
%                     Lcdst = 10*log10(abs(besselh(0,1,k*dst))/abs(besselh(0,1,k*(rri+wi+rs))));
                    Lcdst2 = 20*log10(dst/(rri+wi+rs));
                    appLerch = L_approx_reflbyInput(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);                   
                    app1 = L_approx_reflbyInput1(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);  
                    app2 = L_approx_reflbyInput2(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);  
                    app3 = L_approx_reflbyInput3(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);  
                    app4 = L_approx_reflbyInput4(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);  
                    app5 = L_approx_reflbyInput5(obj,obj.alfa, ws, wi, rs, rri, h1, phiri, lami);  
                    % line source should be Lcdst. But the approximation of the reflection is based on point source. 
                    % Then the correction becomes the point source correction, which should be added on the source. 
                    % Here just put the correction up, there is no physical
                    % meaning
                    approLerch(n, m) = 10*log10(10^(0.1*(appLerch))+10^(0.1*(LbarMff)))+Lcdst2;
                    appro1(n, m) = 10*log10(10^(0.1*(app1))+10^(0.1*(LbarMff)))+Lcdst2;                     
                    appro2(n, m) = 10*log10(10^(0.1*(app2))+10^(0.1*(LbarMff)))+Lcdst2; 
                    appro3(n, m) = 10*log10(10^(0.1*(app3))+10^(0.1*(LbarMff)))+Lcdst2; 
                    appro4(n, m) = 10*log10(10^(0.1*(app4))+10^(0.1*(LbarMff)))+Lcdst2; 
                    appro5(n, m) = 10*log10(10^(0.1*(app5))+10^(0.1*(LbarMff)))+Lcdst2; 
                end
            end
            % the following data is based on point source by timing 1/sqrt(ct)
            % to the line source.
            load([dire  '\' 'levelminusFFL2P.mat'], 'Lrf_oct');
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-approLerch) %#ok<NODEF>
            xlabel('receiver position'); ylabel('frequency'); zlabel('FDTD-approLerch')
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-appro1) % % start from 63Hz
            xlabel('receiver position'); ylabel('frequency'); zlabel('FDTD-appro1')            
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-appro2)
            xlabel('receiver position'); ylabel('frequency'); zlabel('FDTD-appro2')
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-appro3)
            xlabel('receiver position'); ylabel('frequency'); zlabel('FDTD-appro3')
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-appro4)
            xlabel('receiver position'); ylabel('frequency'); zlabel('FDTD-appro4')
            figure()
            surf(rposen(:, 1), frs, Lrf_oct(2:end, :)-appro5)
            xlabel('receiver position'); ylabel('frequency'); zlabel('FDTD-appro5')
        end    
        
        function call_plot(obj)
            close all;
            obj.plot_50imgsrc();
%             obj.plot_approx_50imgsrc();
            obj.plot_1stimag_2to();
%             obj.plot_cmp_FDTD();
%             obj.plot_cmp_FDTD_lerch();
        end
    end
end

