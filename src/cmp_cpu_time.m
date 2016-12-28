classdef cmp_cpu_time
    %compare the running time of single and double diffraction
    
    
    properties
        Xplus = 1.23;
        Xminus = 0.8;
        N = 10000;
    end
    
    methods
        function obj = cmp_cpu_time()
            [tsgTheo, tsgSimpd] = single_diffr_time(obj);
            [tdbTheo, tdbSimpd] = double_diffr_time(obj);
            bar([1,2], [tsgTheo/tsgSimpd, tdbTheo/tdbSimpd], 0.2)
            xlabel('Single or double diffraction')
            ylabel('Time_t_h_e_o_r_y / Time_s_i_m_p_l_i_f_i_e_d')
        end
        
        function [tsgTheo, tsgSimpd] = single_diffr_time(obj)
            tic
            for n=1:obj.N
                DTheo = AD(obj.Xplus)+AD(obj.Xminus);
            end
            tsgTheo = toc;
            tic
            for n=1:obj.N
                DSimpd = 0.37/(0.37+obj.Xplus)+0.37/(0.37+obj.Xminus);
            end
            tsgSimpd = toc;
        end
        
        function [tdbTheo, tdbSimpd] = double_diffr_time(obj)
            tic
            for n=1:obj.N
                DTheo = (AD(obj.Xplus)+AD(obj.Xminus))*(AD(obj.Xplus)+AD(obj.Xminus));
            end
            tdbTheo = toc;
            tic
            for n=1:obj.N
                DSimpd = 0.37/(0.37+obj.Xplus)*0.37/(0.37+obj.Xminus);
            end
            tdbSimpd = toc;
        end
    end
    
end

function vf = f(x)
    [c, s] = fcs(x);
    vf  = (0.5-s)*cos(0.5*pi*x^2)-(0.5-c)*sin(0.5*pi*x^2);
end

function vg = g(x)
    [c, s] = fcs(x);
    vg = (0.5-c)*cos(0.5*pi*x^2)+(0.5-s)*sin(0.5*pi*x^2);
end

function vAD = AD(X)
    vAD = sign(X)*(f(abs(X))-1i*g(abs(X)));
end