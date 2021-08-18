% Fluctuation Correlation Analysis
% Finding txn sites ON, OFF from HybTrack results
% Auto-correlation function calculation
% Larson model fitting

%%%%%%%%%%%%%%%%%%input multiple HybTrack results
%MEF_serum_only_hybtrack_results_Actb or Arc folder
[file,path] = uigetfile('MultiSelect','on');
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
time_int = 8;
N_cell = size(file,2);
s=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
start = [0.1, 20/stemloopsize];%%initial value for Larson fitting
BGThr = 1.5; %%Back ground threshold
fit_res_FCA = [];
on_duration = [];
off_duration = [];
I_TS_on = [];
%%%%%%%%%%%%%%%%%%%%%%%%%ON-OFF detection
for N = 1 : N_cell
    PathFileName=strcat(path,file);
    load(char(PathFileName(1,N)))
    filename = char(file(N));
    filename(end-3:end) = [];
    zerodata = data.tr==0;
    data.tr(zerodata) = NaN;
    N_TS = size(data.tr,2)/3;
    for n = 1 : N_TS
        %%
        I_TS_n_XY=data.tr(:,3*n-2:3*n);
        I_TS_n_XY(isnan(I_TS_n_XY(:,3)),:)=[];
        I_TS_n=I_TS_n_XY(:,3);
        data_size = size(I_TS_n,1);
        BG = mean(I_TS_n_XY(isnan(I_TS_n_XY(:,1)),3));
        if isnan(BG)
            NaN_data = isnan(data.tr);
            NaN_data = [NaN_data(:,size(NaN_data,2)),NaN_data];
            NaN_data(:,size(NaN_data,2)) = [];
            BG = nanmean(data.tr(NaN_data));
        end
        if data_size > 0
            TS_onoff=[NaN;I_TS_n>BGThr*BG;NaN]; 

            on_off_start = [];
            on_off_end = [];
            for i = 1 : size(I_TS_n,1)+1
                if i == 1
                    if TS_onoff(i+1) == 1
                        on_off_start = [on_off_start;[i+1,1]];
                    elseif TS_onoff(i+1) == 0
                        on_off_start = [on_off_start;[i+1,0]];
                    end
                elseif i == size(I_TS_n,1)+1
                    if TS_onoff(i) == 1
                        on_off_end = [on_off_end;[i,1]];
                    else
                        on_off_end = [on_off_end;[i,0]];
                    end
                else
                    if TS_onoff(i) ~= TS_onoff(i+1)
                        if TS_onoff(i+1)==1
                            on_off_start = [on_off_start;[i+1,1]];
                            on_off_end = [on_off_end;[i,0]];
                        elseif TS_onoff(i+1)==0
                            on_off_start = [on_off_start;[i+1,0]];
                            on_off_end = [on_off_end;[i,1]];
                        end
                    end
                end
            end
            on_off_duration = [on_off_end(:,1)-on_off_start(:,1)+1,on_off_start(:,2)];
            on_off_duration(and(on_off_duration(:,2) == 0,on_off_duration(:,1)<8),2) = 1; %short off-times become on-times
            on_off_duration = [[NaN,NaN];on_off_duration;[NaN,NaN]];
            on_off_duration2 = [];

            for i = 1 : size(on_off_duration,1)
                if on_off_duration(i,2) == 1
                    if on_off_duration(i-1,2) ~= 1;
                        on_off_duration2 = [on_off_duration2;[on_off_duration(i,1),1]];
                    else on_off_duration(i-1,2) == 1;
                        on_off_duration2(size(on_off_duration2,1)) = on_off_duration2(size(on_off_duration2,1)) + on_off_duration(i,1);
                    end
                elseif on_off_duration(i,2) == 0
                    if on_off_duration(i-1,2) ~= 0;
                        on_off_duration2 = [on_off_duration2;[on_off_duration(i,1),0]];
                    else on_off_duration(i-1,2) == 0;
                        on_off_duration2(size(on_off_duration2,1)) = on_off_duration2(size(on_off_duration2,1)) + on_off_duration(i,1);
                    end
                end

            end

            on_off_duration2(and(on_off_duration2(:,2) == 1,on_off_duration2(:,1)<8),2) = 0; %short on-times become off-times
            on_off_duration2 = [[NaN,NaN];on_off_duration2;[NaN,NaN]];

            on_off_duration3 = [];
            for i = 1 : size(on_off_duration2,1)
                if on_off_duration2(i,2) == 1
                    if on_off_duration2(i-1,2) ~= 1;
                        on_off_duration3 = [on_off_duration3;[on_off_duration2(i,1),1]];
                    else on_off_duration2(i-1,2) == 1;
                        on_off_duration3(size(on_off_duration3,1)) = on_off_duration3(size(on_off_duration3,1)) + on_off_duration2(i,1);
                    end
                elseif on_off_duration2(i,2) == 0
                    if on_off_duration2(i-1,2) ~= 0
                        on_off_duration3 = [on_off_duration3;[on_off_duration2(i,1),0]];
                    else on_off_duration2(i-1,2) == 0;
                        on_off_duration3(size(on_off_duration3,1)) = on_off_duration3(size(on_off_duration3,1)) + on_off_duration2(i,1);
                    end
                end
            end
            on_off_duration3 = [[NaN,NaN];on_off_duration3;[NaN,NaN]];

            for i = 2 : size(on_off_duration3,1)
                on_off_duration3(i,3) = nansum(on_off_duration3(1:i-1,1))+1;
            end
            
            for i = 2 : size(on_off_duration3,1)
                if on_off_duration3(i,2) == 1
                    on_duration = [on_duration; on_off_duration3(i,1)];
                    I_TS = I_TS_n(on_off_duration3(i,3):on_off_duration3(i,3)+on_off_duration3(i,1)-1);
                    maxsize = max(size(I_TS_on,1),size(I_TS,1));
                    I_TS_on = [[I_TS_on; NaN(abs(maxsize-size(I_TS_on,1)),size(I_TS_on,2))],...
                        [I_TS; NaN(abs(maxsize-size(I_TS,1)),1)]];
                elseif on_off_duration3(i,2) == 0 & on_off_duration3(i-1,2) == 1 & on_off_duration3(i+1,2) == 1
                    off_duration = [off_duration; on_off_duration3(i,1)];
                end
            end
            on_off_duration3(isnan(on_off_duration3(:,1)),:) = [];        
        end
    end
end

%%%%%%%%%%%%%%%%Multi-tau auto-correlation and Larson fitting
for i = 1 : size(I_TS_on,2)
    I_TS_single_on = I_TS_on(:,i);
    I_TS_single_on(isnan(I_TS_single_on))=[];
    if  size(I_TS_single_on,1) > mean(on_duration)%%Data selection:ON traces longer than average ON-time
    ACF_mt = Multi_tau_ACF(I_TS_single_on,s,time_int);
    ACF_mt(size(ACF_mt,1),:)=[];
    x = FCA_fitting(start,ACF_mt,N_stemloop,M_poststemloop);
    fit_res_FCA = [fit_res_FCA;[x(1),(N_stemloop+M_poststemloop)./x(2)]];

    end
end