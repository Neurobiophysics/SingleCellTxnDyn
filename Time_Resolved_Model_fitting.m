function y=IIA_fit_expTt(I_TS_tot,x,X1,X2,L,j)

%specify initial values
X0 = x;
options = optimoptions('lsqnonlin','Display','off','TolFun',1e-8,'FinDiffType','central');
lb=[X0(1)*0.001, X0(2)*0.001,X0(3)*0.001];
ub=[X0(1)*1000, X0(2)*1000, X0(3)*1000];

k_f=lsqnonlin(@(parm)ICS_fun_res(parm,I_TS_tot,X1,X2,L), X0, lb,ub,options);
y=k_f;
simul = figure,
plot(I_TS_tot(:,1),I_TS_tot(:,2)./y(3),'Marker','.','MarkerEdgeColor',[0,0,1],'MarkerSize',24,'LineStyle','none')
hold on
plot(I_TS_tot(:,1),TM(y,I_TS_tot(:,1),X1,X2,L)./y(3),'r','LineWidth',3)
hold on
%text(100, max(TM(y,I_TS_tot(:,1),X1,X2,L))/2,string([string(y(1)),", ", string(y(2))]))
ylabel('normalized amplitude')
xlabel('time (s)')
legend('normalized txn intensity','fitting function')
set(gca,'box','off')
set(gca,'XAxisLocation','bottom',...
'YAxisLocation','left')
set(simul,'WindowStyle','normal')
set(simul,'Position',[-1,-1,30.11*35,30.11*35])
set(gca,'FontSize',41.447)
set(gca,'LineWidth',3)
saveas(gcf,string(j)+'exp'+'.jpg')
saveas(gcf,string(j)+'exp'+'.fig')

%%%%%%%%%%%%%%%%%%%subfunction that computes residuals%%%%%%%%%%%%%%%%%%%%%
function F=ICS_fun_res(x,I_TS_tot,X1,X2,L)

%extract vectors from data, have been passed as parameters
numcol=size(I_TS_tot,2);
tcol=I_TS_tot(:,1);
datacol=I_TS_tot(:,2:1:numcol);
longxcol=repmat(tcol,size(datacol,2),1);
longycol=reshape(datacol,size(datacol,2)*size(datacol,1),1);
F=(TM(x,longxcol,X1,X2,L)-longycol);%./longsigcol;

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