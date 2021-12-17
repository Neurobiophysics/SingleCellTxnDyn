function y=FCA_fitnanmean(x,ACF_mt,N,M)

% Auto-correlation function fitting in the fluctuation correlation analysis
% Model function is based on Larson et al. Science, 2011. 332(6028): p. 475-478.

%specify initial values
X0=x;

%opt=optimset('MaxFunEvals',100000,'MaxIter',5000,'TolFun',1e-10,'TolX',1e-10);
options = optimoptions('lsqnonlin','Display','off');
lb=[X0(1)*0.001, X0(2)*0.001];
ub=[X0(1)*1000, X0(2)*1000];



%[x,resnorm]=lsqnonlin(@(x)ICS_fun_res(x,xmulty,N,M,w), X0, lb,ub,options);
y=lsqnonlin(@(x)ICS_fun_res(x,ACF_mt,N,M), X0, lb,ub,options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subfunction that computes residuals
function F=ICS_fun_res(x,ACF_mt,N,M)

%extract vectors from data, have been passed as parameters
xcol=ACF_mt(:,1);
ycol=ACF_mt(:,2);
errcol=ACF_mt(:,3);
%compute weighted residuals
F=(G(x,xcol,N,M)-ycol)./errcol;
%F=(ICS_fun(x,longxcol)-longycol);%./longsigcol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = G(x,xcol,N,M)
c=x(1);
k=x(2);
G_deno = ((c/k)*(N*(N+1)/2+N*M))^2;
%%%%%%%%%
G_num = zeros(size(xcol));

for n = 1 : N+M
    for m = 1 : n
        log_fact_nm = 0;
        for i = 1 : n-m
            log_fact_nm = log_fact_nm + log(i);
        end
        
        G_num = G_num + I(m,N)*I(n,N)*((c/k))*(exp((n-m).*log(k*xcol)-log_fact_nm-k*xcol));
    end
end
F = G_num/G_deno;



function F = I(n,N)
if n > N
    F = N;
else F = n;
end