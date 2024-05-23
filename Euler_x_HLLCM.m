function [F_ro,F_rou,F_rov,F_E]=Euler_x_HLLCM(ro, rou, rov, E, gamma, k, h)
[f_ro,f_rou,f_rov,f_E,p]=Euler_x(ro,rou,rov,E,gamma);
u=rou./ro;
a=sqrt(gamma*p./ro);
%% Step1, Toro p.332
p_s=1/2*(p(:,1:end-1)+p(:,2:end))-1/2*(u(:,2:end)-u(:,1:end-1)).*...
    (ro(:,1:end-1)+ro(:,2:end)).*(a(:,1:end-1)+a(:,2:end));
p_s(p_s<0)=0;
%% Step2
q_L=sqrt(...
    1+(gamma+1)/(2*gamma)*(p_s./p(:,1:end-1)-1)...
    );
q_L(p_s<=p(:,1:end-1))=1;
S_L=u(:,1:end-1)-a(:,1:end-1).*q_L;
clear("q_L");
q_R=sqrt(...
    1+(gamma+1)/(2*gamma)*(p_s./p(:,2:end)-1)...
    );
q_R(p_s<=p(:,2:end))=1;
S_R=u(:,2:end)+a(:,2:end).*q_R;
clear("q_R","a","p_s");
S=(p(:,2:end)-p(:,1:end-1)+ro(:,1:end-1).*u(:,1:end-1).*(S_L-u(:,1:end-1))-ro(:,2:end).*u(:,2:end).*(S_R-u(:,2:end)))./...
    (ro(:,1:end-1).*(S_L-u(:,1:end-1))-ro(:,2:end).*(S_R-u(:,2:end)));
%% Step3
P_LR=1/2*(p(:,2:end)+p(:,1:end-1)+ro(:,1:end-1).*(S_L-u(:,1:end-1)).*(S-u(:,1:end-1))+...
    ro(:,2:end).*(S_R-u(:,2:end)).*(S-u(:,2:end)));
clear("p");
SL_less=S_L<0;
S_less=S<0;
SR_less=S_R<0;
id_sL=SL_less&(~S_less);
id_sR=S_less&(~SR_less);
id_R=SR_less;
% All fluxes options.
F_ro_sL=S.*(S_L.*ro(:,1:end-1)-f_ro(:,1:end-1))./(S_L-S);
F_rou_sL=(S.*(S_L.*rou(:,1:end-1)-f_rou(:,1:end-1))+S_L.*P_LR)./(S_L-S);
F_rov_sL=S.*(S_L.*rov(:,1:end-1)-f_rov(:,1:end-1))./(S_L-S);
F_E_sL=(S.*(S_L.*E(:,1:end-1)-f_E(:,1:end-1))+S_L.*P_LR.*S)./(S_L-S);
clear("S_L");
F_ro_sR=S.*(S_R.*ro(:,2:end)-f_ro(:,2:end))./(S_R-S);
F_rou_sR=(S.*(S_R.*rou(:,2:end)-f_rou(:,2:end))+S_R.*P_LR)./(S_R-S);
F_rov_sR=S.*(S_R.*rov(:,2:end)-f_rov(:,2:end))./(S_R-S);
F_E_sR=(S.*(S_R.*E(:,2:end)-f_E(:,2:end))+S_R.*P_LR.*S)./(S_R-S);
clear("S_R","S","P_LR");
F_ro=f_ro(:,1:end-1);
F_rou=f_rou(:,1:end-1);
F_rov=f_rov(:,1:end-1);
F_E=f_E(:,1:end-1);
F_ro_R=f_ro(:,2:end);
F_rou_R=f_rou(:,2:end);
F_rov_R=f_rov(:,2:end);
F_E_R=f_E(:,2:end);
%% Choose between fluxes.
F_ro(id_sL)=F_ro_sL(id_sL);
F_rou(id_sL)=F_rou_sL(id_sL);
F_rov(id_sL)=F_rov_sL(id_sL);
F_E(id_sL)=F_E_sL(id_sL);

F_ro(id_sR)=F_ro_sR(id_sR);
F_rou(id_sR)=F_rou_sR(id_sR);
F_rov(id_sR)=F_rov_sR(id_sR);
F_E(id_sR)=F_E_sR(id_sR);

F_ro(id_R)=F_ro_R(id_R);
F_rou(id_R)=F_rou_R(id_R);
F_rov(id_R)=F_rov_R(id_R);
F_E(id_R)=F_E_R(id_R);
end