function [levelMinusFFatL, distanceL, vIdx] = multi_diffr2_theo(spos, rpos, vertex, beta, properClaim, fr, M, PLOT, varargin)
% similar as multi_diffr() but by inputting the source position, receiver
% position and the vertix of the obstacles. 
%
% [levelMinusFFatL, distanceL] = MULTI_DIFFR2(spos, rpos, vertex, beta, properClaim, fr, M)
%   spos is source position as (x, y)
%   rpos is receiver position as (x, y)
%   vertex are the vertex of the obstacles[x1, y1; x2, y2; x3, y3;...; xn, yn]
%   beta is the diffracted angle [beta1, beta2, beta3, ..., betan]
%   properClaim is the piont belongs to a block or to a knife edge shape,
%   as [1, 1, 0, 0, ...] where 1 stands for block and 0 santds for knife
%   edge barriers. properClaim must have the same length with beta.
%   fr is the calculated frequency which can be a number or an array
%   note: the row number of vertex must be the same as the length of beta. 
%   size(vertex, 1)==length(beta)
%   M is the number of adjacent edges
%   PLOT is the switch to plot the path out or not. 1 is switch on.
%
%   levelMinusFFatL is the level minust free field level at distance L
%   distanceL is the reference distance
%
%   [levelMinusFFatL, distanceL]  = MULTI_DIFFR2(spos, rpos, vertix, beta)
%

if size(vertex, 1)~=length(beta)
    error('length of beta must be the same as the row of vertex')
else
    vIdx = shortestPathIndex(spos, rpos, vertex,properClaim, beta);
end
nupi = beta(vIdx==1);
vx = vertex(:,1);
vy = vertex(:,2);
Ldvtx = [spos(1, 1:2); [vx(vIdx==1), vy(vIdx==1)]; rpos(1, 1:2)];
% plot out
if PLOT==1
    figure()
    for n=1:size(vertex, 1)-1
        plot([vertex(n,1) vertex(n,1) vertex(n+1, 1) vertex(n+1, 1)], [0 vertex(n,2) vertex(n+1, 2) 0])
        hold on;
    end
    plot(spos(1), spos(2), 'o', rpos(1), rpos(2), 'o')
    plot(Ldvtx(:, 1), Ldvtx(:, 2), 'k--')
end

propaPath = zeros(size(Ldvtx, 1)-1, 1);
for n=1:size(Ldvtx, 1)-1
    propaPath(n) = dist2D(Ldvtx(n+1,:), Ldvtx(n,:));
end

[thetasn, thetan] = get_angles(spos, rpos, vertex, beta, properClaim, vIdx);
% disp 'path length(in meter): '
% propaPath 
% disp 'thetasn(in degree): '
% thetasn/2/pi*360 
% disp thetan(in degree):
% thetan/2/pi*360 %#ok<*NOPRT>

if nargin==9
    N2xsq = varargin{1};
    [levelMinusFFatL, distanceL] = multi_diffr_theo2(propaPath, nupi, thetasn, thetan, fr, M, N2xsq);
else
    [levelMinusFFatL, distanceL] = multi_diffr_theo2(propaPath, nupi, thetasn, thetan, fr, M);
end

end

function [thetasn, thetan] = get_angles(spos, rpos, vertex, beta, properClaim, edgeIndex)
    vlabl = edge_L_R(beta); % 1 is left edge, 0 is right edge
    thetasn = zeros(sum(edgeIndex), 1);
    thetan = zeros(sum(edgeIndex), 1);
    idx = find(edgeIndex==1);
    for n=1:length(idx)
        if properClaim(idx(n)) == 1 &&  vlabl(idx(n))==1  % left edge of the block
            if n==1 
                thetasn(n) = atan((vertex(idx(n),1)-spos(1,1))/(vertex(idx(n),2)-spos(1,2)));
            else
                thetasn(n) = atan((vertex(idx(n),1)-vertex(idx(n-1),1))/(vertex(idx(n),2)-vertex(idx(n-1),2)));
            end
            if n<length(idx)
                if vertex(idx(n+1),2)-vertex(idx(n),2)>0                    
                    thetan(n) = beta(idx(n))-atan((vertex(idx(n+1),2)-vertex(idx(n),2))/(vertex(idx(n+1),1)-vertex(idx(n),1)));
                else
                    thetan(n) = beta(idx(n));
                end
            else
                thetan(n) = pi+atan((rpos(1,2)-vertex(idx(n),2))/(rpos(11)-vertex(idx(n),1)));
            end
        elseif properClaim(idx(n)) == 1 &&  vlabl(idx(n))==0  % right edge of the block
            if vertex(idx(n),2)-vertex(idx(n-1),2)>0
                error('height of n must be lower than n-1')
            else    
                thetasn(n) = atan(abs(vertex(idx(n),2)-vertex(idx(n-1),2))/(vertex(idx(n),1)-vertex(idx(n-1),1)));
            end
            if n<length(idx)
                thetan(n) = pi+atan(abs(vertex(idx(n+1),2)-vertex(idx(n),2))/(vertex(idx(n+1),1)-vertex(idx(n),1)));
            else
                thetan(n) = pi+atan(abs(rpos(1,2)-vertex(idx(n),2))/(rpos(1,1)-vertex(idx(n),1)));
            end
        else % for knife edge
            if n==1 
                if vertex(idx(n),2)-spos(1,2)<0
                    error('height of n must be higher than source')
                else
                    thetasn(n) = atan((vertex(idx(n),1)-spos(1,1))/(vertex(idx(n),2)-spos(1,2)))-0.5*(2*pi-beta(idx(n)));
                end
            else
                if vertex(idx(n),2)-vertex(idx(n-1),2)<0
                    thetasn(n) = atan(abs(vertex(idx(n),2)-vertex(idx(n-1),2))/(vertex(idx(n),1)-vertex(idx(n-1),1)))+(0.5*pi-0.5*(2*pi-beta(idx(n))));
                else
                    thetasn(n) = atan(abs(vertex(idx(n),1)-vertex(idx(n-1),1))/(vertex(idx(n),2)-vertex(idx(n-1),2)))-0.5*(2*pi-beta(idx(n)));
                end
            end
            if n<length(idx)
                if vertex(idx(n+1),2)-vertex(idx(n),2)<0
                    thetan(n) = 0.5*pi-0.5*(2*pi-beta(idx(n)))+pi+atan(abs(vertex(idx(n+1),2)-vertex(idx(n),2))/(vertex(idx(n+1),1)-vertex(idx(n),1)));
                else
                    error('height of n+1 must be lower than n')                    
                end         
            else
                thetan(n) = 0.5*pi-0.5*(2*pi-beta(idx(n)))+pi+atan(abs(rpos(1,2)-vertex(idx(n),2))/(rpos(1,1)-vertex(idx(n),1)));
            end
        end
    end    
end

function vIdx = shortestPathIndex(spos, rpos, vertex,properClaim, beta)
%  to get the edges indexes of the shortest diffraction path 
%   VALID FOR 2.5D ONLY, i.e. FLAT ROOF AND KNIFE EDGE BARRIERS
%   suppose the building and barriers are on the ground, i.e. the y==0 for
%   all the buildings and blocks and the terrain effect  is NOT included. 
%
% vIdx = shortestDiffrPath(spos, rpos, vertex)
%   spos is source position as (x, y)
%   rpos is receiver position as (x, y)
%   vertex are the vertex of the obstacles[x1, y1; x2, y2; x3, y3;...; xn, yn]
%   properClaim is the piont belongs to a block or to a knife edge shape,
%   as [1, 1, 0, 0, ...] where 1 stands for block and 0 santds for knife
%   edge barriers. properClaim must have the same length with the row of vertex.
%   beta is the diffracted angle [beta1, beta2, beta3, ..., betan]
%
%   vIdx is the index of the vertex, as vIdx = [1 0 1 1 0 1]

gsrv = [spos(1,1:2); vertex; rpos(1, 1:2)];
Ldvtx = zeros(size(vertex,1), 2);
LR = edge_L_R(beta);
p = 1;
cnt =1;
while p<=size(gsrv,1)-2
    M = gsrv(p,:);
    MptN = gsrv(p+1, :);
    q = p+2;
    Q = q;
    while q<=size(gsrv,1)
        MptN2 = gsrv(q-1, :);
        N = gsrv(q, :);            
        askInter = polyxpoly([M(1), N(1)], [M(2), N(2)], [MptN(1), MptN(1)], [MptN(2), 0]);            
        askInter2 = polyxpoly([M(1), N(1)], [M(2), N(2)], [MptN2(1), MptN2(1)], [MptN2(2), 0]);
        if q-p>=4
            for t=(p+2):(q-2)   % p+1 is MptN, q-1 is MptN2         
                MptN3 = gsrv(t, :);
                askInter3 = polyxpoly([M(1), N(1)], [M(2), N(2)], [MptN3(1), MptN3(1)], [MptN3(2), 0]);
                if isempty(askInter)==1 && isempty(askInter2)==1 && isempty(askInter3)==1
                    Q = q;
                end
            end
        else
            if isempty(askInter)==1 && isempty(askInter2)==1     
                Q = q;
            end 
        end        
        q = q+1; 
    end 
    if Q<=size(vertex,1)+1  % Q is not the receiver
        if properClaim(Q-2)==0 || LR(Q-2)==1  % left edge or thin barrier
            if Q==3 % p is the source position
                Ldvtx(cnt, :) = gsrv(2,:);
                p = 2;
            else
                Ldvtx(cnt, :) = gsrv(Q,:);
                p = Q;
            end
        else % right edge of a block
            askInter3 = polyxpoly([M(1), gsrv(Q,1)], [M(2), gsrv(Q,2)], [gsrv(Q-1,1), gsrv(Q-1,1)], [gsrv(Q-1,2), 0]); 
            if isempty(askInter3)==1
                Ldvtx(cnt, :) = gsrv(Q,:);
                p = Q;
            else
                Ldvtx(cnt, :) = gsrv(Q-1,:);
                p = Q-1;
            end
        end        
    else % Q is the receiver
        askInter4 = polyxpoly([Ldvtx(cnt-1,1), gsrv(Q,1)], [Ldvtx(cnt-1,2), gsrv(Q,2)], [gsrv(Q-1,1), gsrv(Q-1,1)], [gsrv(Q-1,2), 0]); 
        if isempty(askInter4)==0
            Ldvtx(cnt, :) = gsrv(Q-1,:);
        end
        p = Q-1;
    end
    cnt = cnt+1;  
end
rLdvtx = Ldvtx(:,1);
idx = find(rLdvtx==0);
idx = idx(1);
if idx~=1
    Ldvtx2 = Ldvtx(1:idx-1,:);
else
    error('Only one point in the diffracted edge variable!!!')
end
vIdx = zeros(size(vertex, 1), 1);
for n=1:size(Ldvtx2,1)
    idx = vertex(:, 1)==Ldvtx2(n,1);
    vIdx(idx) = 1;
end
end

function LR = edge_L_R(beta)
% get the edge is left or right edge of a block
% 
% LR = edge_L_R(beta)
%   beta is the diffraction angle [1.5pi, 1.5pi, ....]
%   LR is the sign of left or right edge,[1 0 1 0 1 0 99 1 ...]
%   1 is left edge, 0 is right edge, 99 is a thin barrier

    edpos = find(beta==2*pi);    
    if isempty(edpos)==1
        LR = rem(1:length(beta),2);
    else
        labels = zeros(length(beta)-length(edpos), 1);
        n=1;
        m=1;
        while n<=length(beta)
            g = isempty(find(edpos==n, 1));
            if g==1  % empty
                labels(n) = m;
            else
                labels(n) = 99; %99 stands for edge
                m = 0;
            end
            n = n+1;
            m = m+1;
        end
        LR = zeros(length(beta),1);
        for p=1:length(beta)
            if labels(p)~=99
                LR(p) = rem(labels(p), 2); 
            else
                LR(p) = 99;
            end
        end
    end 
end

function d = dist2D(p1, p2)
    d = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
end
    
    
    
    