% Estimation of the GFE model with random coefficients (coefficient of lagged income is heterogeneous)

clear all;
clc;

load('final_data');

X0 = ldem_linc;
Y0 = dem;

XX = X0(:,1:2);
YY = Y0;
clear X Y
X = XX;
Y = YY;
T = 7;
K = 2;
N = size(Y,1)/T;

% Deviations to the mean (optional)
%Ybar=zeros(N,1);
%Xbar=zeros(N,K);
%for i=1:N
%    Ybar(i)=mean(Y((i-1)*T+1:i*T));
%    Xbar(i,:)=mean(X((i-1)*T+1:i*T,:));
%end
   
%Y=Y-kron(Ybar,ones(T,1));
%X=X-kron(Xbar,ones(T,1));

G = 3;

gi_auxaux=zeros(N,G);
gi=zeros(N,G);
gitot=zeros(N*T,G*T);
Xgi=zeros(N*T,G*K);
Nsim=1000;
Resbeta=zeros(Nsim,G*K);
Resdelta=zeros(Nsim,G*T);
ResQ=zeros(Nsim,1);


for jsim=1:Nsim

    % Starting values: g_i
    gi_init=zeros(N,G);
    for i=1:N
        perm=randperm(G);
        for g=1:G
            gi_init(i,g)=(perm(1)==g);
        end
    end
    
    gi=gi_init;
    
    par_init=zeros(G*(K+T),1);
    
    % iterations
    deltapar=1;
    while deltapar>0
       
        % step 1: update
     
        for i = 1:N
            giaux=kron(eye(T),gi(i,:));
            gitot((i-1)*T+1:T*i,:) = giaux;
            Xgi((i-1)*T+1:T*i,:) = kron(X((i-1)*T+1:i*T,:),gi(i,:));
        end
        
        Xtot=[gitot Xgi];
        par_new=pinv(Xtot)*Y;
        
        delta_new_vect=par_new(1:G*T);
        delta_new=reshape(delta_new_vect,G,T)';
        beta_new_vect=par_new(G*T+1:G*T+G*K);
        beta_new=reshape(beta_new_vect,G,K)';
        
        % step 2: assignment
        
        for g = 1:G
            U = (Y-X*beta_new(:,g)-kron(ones(N,1),delta_new(:,g))).^2;
            RU = reshape(U,T,N)';
            gi_auxaux(:,g) = sum(RU')';
        end;
        gi_class = min(gi_auxaux')';
        
        for g = 1:G
            gi(:,g) = (gi_auxaux(:,g) == gi_class);
        end
        
        
        
        deltapar=norm(par_new-par_init);
        
        par_init=par_new;
        
        
        %(Y-Xtot*par_init)'*(Y-Xtot*par_init)
    end
    
    Resbeta(jsim,:)=beta_new_vect';
    Resdelta(jsim,:)=delta_new_vect';
    ResQ(jsim)=(Y-Xtot*par_new)'*(Y-Xtot*par_new);
    jsim
  

end


[x1,x2]=sort(ResQ,'ascend');
%Resbeta(x2,:)
%plot(Resbeta(x2,1),Resbeta(x2,2),'.r')

beta_final=Resbeta(x2(1),:)';
delta_final=Resdelta(x2(1),:)';

par_final=[delta_final;beta_final];

delta_final = reshape(delta_final,G,T)';
beta_final = reshape(beta_final,G,K)';



% assignment

for g = 1:G
    U = (Y-X*beta_final(:,g)-kron(ones(N,1),delta_final(:,g))).^2;
    RU = reshape(U,T,N)';
    gi_auxaux(:,g) = sum(RU')';
end;
gi_class = min(gi_auxaux')';

for g = 1:G
    gi(:,g) = (gi_auxaux(:,g) == gi_class);
end

gi0=gi;

Inc=zeros(N,T);
for t=1:T
    Inc(:,t)=X(t:T:N*T,2);
end

vecgi=gi*(1:1:G)';
Matinc=zeros(G,T);
for g=1:G
    Matinc(g,:)=mean(Inc((vecgi==g),:));
end

plot(Matinc')

Dem=zeros(N,T);
for t=1:T
    Dem(:,t)=Y(t:T:N*T,1);
end

Matdem=zeros(G,T);
for g=1:G
    Matdem(g,:)=mean(Dem((vecgi==g),:));
end

plot(Matdem')

% Point estimates

beta_final
beta_final(2,:)./(1-beta_final(1,:))

% Standard errors (large-N,T)

for i = 1:N
    giaux=kron(eye(T),gi(i,:));
    gitot((i-1)*T+1:T*i,:) = giaux;
    Xgi((i-1)*T+1:T*i,:) = kron(X((i-1)*T+1:i*T,:),gi(i,:));
end

Xtot=[gitot Xgi];

        
ei=Y-Xtot*par_final;
Rei = reshape(ei,T,N)';
Omega = zeros(N*T,N*T);
Mi = zeros(T,T);
for i = 1:N
   Mi = Rei(i,:)'*Rei(i,:);
   Omega((i-1)*T+1:i*T,(i-1)*T+1:i*T) = Mi;
end

V = inv(Xtot'*Xtot)*Xtot'*Omega*Xtot*inv(Xtot'*Xtot);
V=V*N*T/(N*T-T-K);
std_cluster = sqrt(diag(V));

std_cluster(G*T+1:G*T+K*G)

std_theta_inc=zeros(G,1);

for g=1:G
    Mat_inc=[beta_final(2,g)/(1-beta_final(1,g))^2 1/(1-beta_final(1,g))];
    Var_theta_inc=Mat_inc*V([G*T+g,G*T+G+g],[G*T+g,G*T+G+g])*Mat_inc';
    std_theta_inc(g)=sqrt(Var_theta_inc);
end

std_theta_inc
