%% Time-Resolved Measurement Simulation
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
elong_rates = [10^0.5 ,10, 10^1.5,100]; % elongation rate input
T_terms = [1,10,100,250,500]; % termination time input
for ini_rate = ini_rates
        for elong_rate = elong_rates
            for T_term = T_terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Single Simulation Input %%%%%%%%%%%%%%%%
%ini_rate = 0.12; % initiation rate input
%elong_rate = 13.0; % elongation rate input
%T_term = 228; % termination rate input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            [ini_rate, elong_rate, T_term]
I_TS_simul = [];
I_TS_tot = [];
fit_res = [];
for n = 1 : 1000
    time_tot = frm_max * time_int; % total observation time
    N_rough = round((L/elong_rate) * ini_rate)*2;
    txn_ini_time = cumsum(exprnd(1/ini_rate,1,2*N_rough)); % initiation with Poissonian promoter
    bp_0t_elong = poissrnd(elong_rate*((L/elong_rate) - txn_ini_time), size(txn_ini_time));
    bp_0t_elong = bp_0t_elong(bp_0t_elong<L);
    N_elong = size(bp_0t_elong,2);
    %N_elong = round(L*(ini_rate/elong_rate));
    %bp_0t_elong = L*rand(1,N_elong);
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
    I_TS_simul = [I_TS_simul,I_TS];
    I_TS = [[0:time_int:(size(I_TS)-1)*time_int]',I_TS]; 

    
    N1 = size(I_TS(:,2),1);
    t_tot = time_int*N1;

    fit_res_temp = Time_Resolved_Model_fitting(I_TS,[20,100,max(I_TS(:,2))],X1,X2,L);
    if fit_res_temp(1)<1000; %remove data which has no plateau
    fit_res = [fit_res;fit_res_temp];
    I_TS_tot = [I_TS_tot,I_TS(:,2)];
    else
    end
    if size(fit_res,1)>=n_sample
        break
    else
    end
end
tcol = [0:time_int:(499)*time_int]';

%%
%y=mean(fit_res)
%avg_I_TS = mean(I_TS_tot./repmat(fit_res(:,3)',500,1),2);
%err_I_TS = std((I_TS_tot./repmat(fit_res(:,3)',500,1))');
%fit_I_TS = TM(y,tcol,X1,X2,L)./y(3);
%simul = figure,
%hold on
%plot(tcol(1:251),avg_I_TS(1:251),'Marker','.','MarkerEdgeColor',[0,0,1],'MarkerSize',24,'LineStyle','none')
%errorbar(tcol(1:251),avg_I_TS(1:251), err_I_TS(1:251),'b')
%plot(tcol(1:251),fit_I_TS(1:251),'r','LineWidth',3)


%ylabel('normalized amplitude')
%xlabel('time (s)')
%xlim([0 2000])
%ylim([-0.1 1.2])
%set(gca,'box','off')
%set(gca,'XAxisLocation','bottom',...
%'YAxisLocation','left')
%set(simul,'WindowStyle','normal')
%set(simul,'Position',[-1,-1,30.11*35,30.11*35])
%set(gca,'FontSize',41.447)
%set(gca,'LineWidth',3)
%saveas(gcf,string(j)+'exp'+'.jpg')
%saveas(gcf,string(j)+'exp'+'.fig')


fit_res_mean = mean(fit_res);
fit_res_std = std(fit_res);
fit_res_mean(3) = [];
ci_ki_Tti_kr_Ttr = [ci_ki_Tti_kr_Ttr;[ini_rate,elong_rate,T_term,fit_res_mean]];
%%%%%%%%%%%%%%%%%% Multiple Simulation For-loop end%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Model function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y1 = TM(x,tcol,X1,X2,L)

k = x(1);
T_term = x(2);
I_M = x(3);
t_i = 0;
y1 = zeros(size(tcol));
for i = 1 : size(tcol,1)
    if i == 1
        y1(i) = I_M;


    elseif tcol(i)-t_i-L/k < 0
        y1(i) = I_M*(heaviside((tcol(i)-t_i))*heaviside(X1/k-(tcol(i)-t_i)) ...
            + heaviside((tcol(i)-t_i)-X1/k)*heaviside(X2/k-(tcol(i)-t_i))*(1-((k*(tcol(i)-t_i)-X1)^2)/((2*(L+k*T_term)-X1-X2)*(X2-X1))) ...
            + heaviside((tcol(i)-t_i)-X2/k)*heaviside(L/k-(tcol(i)-t_i))*(L+k*T_term-k*(tcol(i)-t_i))/(L+k*T_term-(X1+X2)/2));
    else
        y1(i) =  I_M*(heaviside((tcol(i)-t_i)-L/k)*((2*k*T_term*exp(-((tcol(i)-t_i-L/k)/T_term)))/(2*(L+k*T_term)-X1-X2)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end