%%%%%%%%%%%%%%%%%% Load gene structure input
%gene_structure_Actb or Arc.mat
[genefile,genepath] = uigetfile('*.mat');
genePathFileName=strcat(genepath,genefile);
load(char(genePathFileName));
%%%%%%%%%%%%%%%%%% Gene Sturucture Input
N_stemloop =		input_data.N; % length of stem-loop region (number of stem-loops)
M_poststemloop =	input_data.M; % length of post-stem-loop region (scaled by stem-loop size)
X1          = 		input_data.X1; % 5' end position of stem-loop region (unit: bps)
X2          = 		input_data.X2; % 3' end position of stem-loop region (unit: bps)
L           = 		input_data.L; % 3' end position of the gene (unit: bps)
stemloopsize =  	input_data.stemloopsize; % stem-loop size (unit: bps)
stemloopend =   	input_data.stemloopend; % 3'end position of each stem-loop (unit: bps)
%%%%%%%%%%%%%%%%%%
frm_max = 500; % The number of total measurement frames
time_int = 8; % Time interval of microscope measurement
n_sample = 100; % number of samples
%%
ci_ki_Tti_kr_Ttr = []; % c_input, k_input, Tt_input, k_output, Tt_output
%%%%%%%%%%%% Multiple Simulation Input %%%%%%%%%%%%%%%%
ini_rates = [0.01,10^(-1.5),0.1,10^(0.5),1]; % initiation rate input 
elong_rates = [10^0.5 ,10, 10^1.5,100, 10^2.5]; % elongation rate input
T_terms = [1,10,100,250,500]; % termination time input
for ini_rate = ini_rates
        for elong_rate = elong_rates
            for T_term = T_terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Single Simulation Input %%%%%%%%%%%%%%%%
%ini_rate = 0.1; % initiation rate input
%elong_rate = 10; % elongation rate input
%T_term = 250; % termination rate input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            [ini_rate, elong_rate, T_term]

fit_res = [];
for n = 1 : 1000
    N_elong = round(L*(ini_rate/elong_rate));
    bp_0t_elong = L*rand(1,N_elong);
    Tt_elong = exprnd(T_term,1,N_elong);
    bp_t_elong = cumsum([bp_0t_elong; poissrnd(elong_rate*time_int,frm_max-1,N_elong)]);
    bp_max_elong = elong_rate*Tt_elong+L;
    bp_t_elong(bp_t_elong>repmat(bp_max_elong,size(bp_t_elong,1),1)) = NaN;
    I_TS_cumsum_elong = cumsum(bp_t_elong>reshape(stemloopend,1,1,N_stemloop),3);
    I_TS_elong = I_TS_cumsum_elong(:,:,N_stemloop);
    
    N_term = round(ini_rate*T_term);
	bp_0t_term = exprnd(T_term,1,N_term);
    bp_t_term = cumsum([bp_0t_term;(-time_int)*ones(frm_max-1,N_term)]);
    bp_t_term(bp_t_term>0) = N_stemloop;
    bp_t_term(bp_t_term<=0) = 0;
    I_TS = sum([I_TS_elong,bp_t_term],2,'omitnan');
    I_TS = I_TS + normrnd(1,1,[size(I_TS,1),1]);
    I_TS = [[time_int:time_int:size(I_TS)*time_int]',I_TS]; 

    
    N1 = size(I_TS(:,2),1);
    t_tot = time_int*N1;

    fit_res_temp = IIA_fitting(I_TS,[20,100,max(I_TS(:,2))],X1,X2,L);
    if fit_res_temp(1)<1000; %remove data which has no plateau
    fit_res = [fit_res;fit_res_temp];
    else
    end
    if size(fit_res,1)>=n_sample
        break
    else
    end
end

fit_res_mean = mean(fit_res);
fit_res_mean(3) = [];
ci_ki_Tti_kr_Ttr = [ci_ki_Tti_kr_Ttr;[ini_rate,elong_rate,T_term,fit_res_mean]];
%%%%%%%%%%%%%%%%%% Multiple Simulation For-loop end%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
