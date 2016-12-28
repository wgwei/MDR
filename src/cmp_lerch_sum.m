% cmp_lerch_sum
% compare the calculating efficiency between Lerch transcendent and
% summation

% prepare inputs
spos = [0 0 0];
rpos = [30 0 0];
N1p = [10 10];
N2p = [20 10];
fr = 500;
betas = [pi/4];
betar = [pi/4];

Wi = 10;
wavlen = 1.7;
h1 = 10;
rr = 15;
z = 0.97;
s = 2;
Ws = 10;
a = (sqrt(2*Wi/wavlen)*h1+0.5*Ws+rr+Wi+Ws)/Ws;

N = 100;
M = [5 10 20 30 40 50 60 70 80 90 100];

avg_summation = zeros(100, 1);
avg_summation2 = zeros(100, 1);
avg_lerchv = zeros(100, 1);
ratio_plot = zeros(length(M), 1);

% run comparison
for p=1:length(M)
    MM = M(p);
    for r=1:N
        tic
        for n=1:MM
            [p2lvlSimpd, p2lvlTheo, L] = double_diffr_for_cmp_lerch(spos, rpos, N1p, N2p, fr, betas, betar);            
        end
        avg_summation(r) = toc;  

%         tic
%         for n=1:MM
%             [p2lvlSimpd, p2lvlTheo, L] = double_diffr_for_cmp_lerch2(spos, rpos, N1p, N2p, fr, betas, betar);
%         end
%         avg_summation2(r) = toc;
        
        tic
        value = lerch ( z, s, a );
        avg_lerchv(r) = toc;
    end

    % output ratio    
    ratio_plot(p, 1) = mean(avg_summation)/mean(avg_lerchv)
%     ratio_plot(p, 2) = mean(avg_summation2)/mean(avg_lerchv)
end
plot(M,ratio_plot, '-o')
xlabel('Number of image sources')
ylabel('Time_s_u_m / Time_L_e_r_c_h')

