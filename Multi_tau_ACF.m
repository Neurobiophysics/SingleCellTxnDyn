function ACF_mt = MultiTauACF(I_TS,s,time_int)
% Auto-correlation function by multi-tau method based on
% Wohland, T., R. Rigler, and H. Vogel, The standard deviation in fluorescence correlation spectroscopy.
% Biophysical journal, 2001. 80(6): p. 2987-2999.

n_chs1 = 2^s;
frm_max = size(I_TS,1);
p = frm_max;
n = 1;
while p > 1
    n = n + 1;
    p = frm_max/2^n;
end
N_chgroup = n-s+1;
N_ch = (n_chs1/2)*(N_chgroup+1);
chstart = transpose([1,n_chs1+1:n_chs1/2:(N_chgroup+1)*(n_chs1/2)+1]);
ch = zeros(N_ch,frm_max);
ACF_frm = zeros(N_ch,frm_max);



frm_ch = (0:1:2^s-1);
for i = 1 : size(chstart,1)-2
    frm_ch = [frm_ch,2^(s+i-1):2^(i):(2^(s+i)-1)];
end
frm_ch = transpose(frm_ch);

frm_ch_1 = [frm_ch;frm_ch(end)];
frm_ch_1(1) = [];
time_int_ch = frm_ch_1-frm_ch; %Channel time intervals
time_int_ch(end) = time_int_ch(end-1);
max(frm_ch(frm_ch < size(I_TS,1)));
I_TS(max(frm_ch(frm_ch < size(I_TS,1)))+1:size(I_TS,1)) = [];
frm_max = size(I_TS,1);

multitau = (0:time_int:(2^s-1)*time_int);
for i = 1 : size(chstart,1)-2
    multitau = [multitau,2^(s+i-1)*time_int:2^(i)*time_int:(2^(s+i)-1)*time_int];
end
multitau = transpose(multitau);
n_dir = [];
n_del = [];

frm_ch(frm_ch>size(I_TS,1))=[];
for i_tau = 1 : size(frm_ch,1)
    tau = frm_ch(i_tau);
    ch_temp = movsum(I_TS,time_int_ch(i_tau))';
    
    n_dir_temp = ch_temp(floor(time_int_ch(i_tau)/2)+1:end-floor(time_int_ch(i_tau)/2)-tau);
    n_dir_temp = n_dir_temp(1:time_int_ch(i_tau):end);
	n_dir(i_tau,:) = [n_dir_temp, NaN(1,size(I_TS,1)-size(n_dir_temp,2))];%Measured intensities from direct monitor by tau (row) and time (column)
    
    n_del_temp = ch_temp(floor(time_int_ch(i_tau)/2)+1+tau:end-floor(time_int_ch(1)/2));
    n_del_temp = n_del_temp(1:time_int_ch(i_tau):end);
    n_del(i_tau,:) = [n_del_temp, NaN(1,size(I_TS,1)-size(n_del_temp,2))];%Measured intensities from delay monitor by tau (row) and time (column)
end

M = repelem(frm_max,2^s); 
for i = 1 : size(chstart,1)-2
    M = ([M,repelem(fix(frm_max/2^(i)),2^(s-1))]);
end
M = transpose(M);
m = transpose([0:2^s-1,repmat(2^(s-1)+1:1:2^s,1,size(chstart,1)-2)]); % del
M_m = M(1:size(n_dir,1))-m(1:size(n_dir,1));% M-m is the number of possible measurements by the channels

M_dir = sum(n_dir,2,'omitnan')./M_m;%average intensity over time from direct monitor
M_del = sum(n_del,2,'omitnan')./M_m;%average intensity over time from direct monitor

ACF = sum((n_dir-repmat(M_dir,1,size(n_dir,2))).*(n_del-repmat(M_del,1,size(n_del,2))),2,'omitnan')./(M_m.*M_dir.*M_del);
%%%%%%%%Saffarian's SEM Saffarian Biophysical journal, 2003. 84(3): p. 2030-2042.
ACF_SEM...
    = sqrt(...
    sum(((n_dir-repmat(M_dir,1,size(n_dir,2))).^2).*((n_del-repmat(M_del,1,size(n_del,2))).^2),2,'omitnan')./M_m...
    - sum((n_dir-repmat(M_dir,1,size(n_dir,2))).*(n_del-repmat(M_del,1,size(n_del,2))),2,'omitnan').^2./(M_m.^2)...
    )...
./(sqrt(M_m).*M_dir.*M_del);


ACF(isnan(ACF))=[];

ACF_mt = [multitau(1:size(ACF,1))+0.1*ones(size(ACF)),ACF,ACF_SEM(isfinite(ACF_SEM))];
