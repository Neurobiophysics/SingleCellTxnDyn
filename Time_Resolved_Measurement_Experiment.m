%% Time-Resolved Measurement
% Input intensity time traces and gene structure
%%%%%%%%%%%%%%%%%% Input selected HybTrack results
%IIA_selected_TS_Actb or Arc.mat
[file,path] = uigetfile('.mat');
%%%%%%%%%%%%%%%%%% Load selected Hybtrack results 
PathFileName=strcat(path,file);
load(char(PathFileName))
%%%%%%%%%%%%%%%%%% load gene structure input
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
time_int = 8; % Microsope imaging time interval
%%
fit_res_IIA = [];
for ii = 1 : size(I_TS_selected,2)
    I_TS = [I_TS_selected(:,ii)];
    ini_val = [5,100,I_TS(1)]; % Initial value for k, T_t, and F(0)
    tcol = transpose([0:time_int:time_int*(size(I_TS,1)-1)]); % Generating time axis
    I_TS_t = [tcol,I_TS];
    I_TS_t(any(isnan(I_TS_t), 2), :) = [];
    res = Time_Resolved_Model_fitting(I_TS_t,ini_val,X1,X2,L,ii); %exponential termination time
    fit_res_IIA = [fit_res_IIA;res];
end