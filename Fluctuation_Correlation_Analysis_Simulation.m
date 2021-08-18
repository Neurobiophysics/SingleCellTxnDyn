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
s=4; % multi-tau channel order


%%
ci_ki_Tti_cr_Tdr = []; % c_input, k_input, Tt_input, c_output, Td_output

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


c = ini_rate;
k = elong_rate;
ini_val = [c,k/stemloopsize];
fit_res = [];
for n = 1 : 100 % number of samples
    N_txn = round(frm_max * 1 * time_int * ini_rate); % number of mRNA Txn within total time
    bp = NaN(frm_max,N_txn);    
    time_tot = frm_max * time_int; % total observation time
    txn_ini_time = (time_tot*1)*rand(1,N_txn); % times that each mRNA Txn initiated
    txn_ini_frm = ceil(txn_ini_time/time_int);
    bp_txn_ini_frm = poissrnd((ceil(txn_ini_time/time_int)-txn_ini_time/time_int)*elong_rate);
    %txn_term_time = exprnd(T_term,1,N_txn);
    txn_term_time = T_term*ones(1,N_txn);
    bp_txn = cumsum([bp_txn_ini_frm;poissrnd(elong_rate*time_int,frm_max-1,N_txn)]);
    bp_tot_extended = elong_rate*txn_term_time + L;
    repmat(bp_tot_extended,size(bp_txn,1),1);
    bp_tot_extended_temp = repmat(bp_tot_extended,size(bp_txn,1),1);
    bp_txn(bp_txn >= bp_tot_extended_temp) = NaN;

    
    
    
    frm_extend = max(rem(find(~isnan(bp_txn)),frm_max));
    bp_txn(frm_extend+1:size(bp_txn,1),:)=[];
    bp = [bp;zeros(frm_extend,N_txn)];
    bp(reshape(cumsum([txn_ini_frm;ones(frm_extend -1,size(txn_ini_frm,2))]) + repmat([0:size(bp,1):size(bp,1)*(size(bp,2)-1)],frm_extend,1),frm_extend*N_txn,1)) = reshape(bp_txn,frm_extend*N_txn,1);
    bp(frm_max+1:size(bp,1),:) = [];
    bp(bp<stemloopend(1)) = NaN;
    for i = 1 : size(stemloopend,1)
        bp(bp>=stemloopend(size(stemloopend,1)-i+1)) = size(stemloopend,1)-i+1;
    end
 
    I_TS = sum(bp','omitnan')';

    %%%%%%%%%%%
    corr_start = find(I_TS ~= 0, 1, 'first');
    corr_end = find(I_TS~=0,1,'last');
    I_TS = I_TS(corr_start:corr_end);
    %I_TS = I_TS(100:312);
    %%%%%%%%%%%%% Data crop, excluding the zeros from both ends
    I_TS = I_TS + normrnd(1,1,[size(I_TS,1),1]);
    I_TS = [transpose(time_int:time_int:size(I_TS)*time_int),I_TS];    
      
    N1 = size(I_TS(:,2),1);
    t_tot = time_int*N1;
 
    ACF_mt = Multi_tau_ACF(I_TS(:,2),s,time_int); % Calculating auto-correlation function with multitau method
    ACF_mt(any(ACF_mt==0,2),:) = [];
 
    fit_res_temp = FCA_fitting([ini_val(1),ini_val(2)],ACF_mt,N_stemloop,M_poststemloop); %Larson model fitting

    
    fit_res = [fit_res;fit_res_temp];
  
fit_res_mean = mean(fit_res);

%%%%%%%%%%%%%%%%%% Multiple Simulation For-loop end%%%%%%%%%%%%%%%%%%%%%
end
ci_ki_Tti_cr_Tdr = [ci_ki_Tti_cr_Tdr;[ini_rate,elong_rate,T_term,fit_res_mean(1),(N_stemloop+M_poststemloop)/fit_res_mean(2)]];
%%%%%%%%%%%%%%%%%% Multiple Simulation For-loop end%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%