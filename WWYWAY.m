function [Z,t,W]=WWYWAY(X,lambda_c,lambda_d,mu,epsilon,max_time,gamma,Noise_ratios)%Noise_ratios为X中的噪声比
V=numel(X);
N=size(X{1},2);
C=zeros(N,N);
Y22=zeros(N,N);
G=zeros(N,N);
sigma=zeros(N,1);%
Y=cell(V,1);
A=cell(V,1);
W=cell(V,1);
D=cell(V,1);
for i=1:V
    Y{i}=zeros(size(X{i}));
    W{i}=ones(size(X{i}))*(1/(size(X{i},1)*N));
    D{i}=rand(N);
    A{i}=X{i}-X{i}*(C+D{i});
end
t=0;
while t<max_time
    C_L=C;
    G_L=G;
    D_L=D;
    
    t=t+1;
    %% update G 
    AA=C-(Y22/mu);
    [U,S,V2]=svd(AA);
    a=zeros(N,1);
    k=1;
    while(k<=1)%
      for j=1:N
          a(j)=(exp(-(sigma(j)/gamma)*((sigma(j)/gamma)+2)))*((sigma(j)+gamma)/(gamma*gamma))*2;
      end
      sigma=diag(S)-(lambda_c/mu)*a;
      sigma(sigma<0)=0;
      k=k+1;
    end
    G=U*diag(sigma)*V2';
    
    %% update C 
    C_a=eye(N);
    C_b=G+(Y22/mu);
    for j=1:V
        C_a=C_a+X{j}'*X{j};
        C_b=C_b+X{j}'*(X{j}-A{j}-(X{j}*D{j})-(Y{j}/mu));
    end
    C=C_a\C_b;
    
    %% update D
    for i=1:V
        left_d=2*lambda_d*eye(N)+mu*X{i}'*X{i};
        right_d=mu*X{i}'*(X{i}-A{i}-(X{i}*C)-(Y{i}/mu));
        D{i}=left_d\right_d;
    end
    
    %% update A
    for i=1:V
        left_a=W{i}+mu*ones(size(X{i}));
        right_a=mu*X{i}-mu*X{i}*(C+D{i})-Y{i};
        A{i}=right_a./left_a;  
    end
    
     %% update W
    for i=1:V
        ss=size(X{i},1);
        R_T=A{i}.*A{i};
        b=R_T(:);
        B=sort(b);
        R=floor(Noise_ratios(i)*ss*N);
        
        c=0;
        for j=1:(ss*N-R)
            c=c+B(j);
        end
        
        W{i}= ((2*B(ss*N-R+1))-R_T)/((2*B(ss*N-R+1)*(ss*N-R))-c);
        W{i}(W{i}<0)=0;
    end
    
    %% update Y
    for i=1:V
        Y{i}=Y{i}+mu*(A{i}-X{i}+X{i}*(C+D{i}));
    end
    %% update Y22
    Y22=Y22+mu*(G-C);
    %% update mu
    mu_max=10^6;
    rho=1.5;
    mu=min(mu_max,rho*mu);
    
    %% check if converge
    diff1=norm(C_L-C,inf);
    diff2=norm(G_L-G,inf);
    diff3=zeros(V,1);
    for i=1:V
        diff3(i)=norm(D_L{i}-D{i},inf);
    end
    diff3_max=max(diff3);
    if (max([diff3_max diff1 diff2]))<epsilon
       break;
    end
end
Z=(abs(C)+abs(C'))/2;
for i=1:V
    Z=Z+(abs(D{i})+abs(D{i}'))/(2*V);
end
end