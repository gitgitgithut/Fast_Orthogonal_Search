close all, clear all, clc;

% initialization
N = 1000;
fun = '1+0.7x(n)+0.8x(n-1)+0.1x(n-1)y(n-1)';
sigma_x = 1;
mu_x = 0;
x = sigma_x*randn(1,3*N)+mu_x;
y = zeros(1,3*N);
for n= 1:3*N
    y(n) = 1+0.7*delay(x,n,0)+0.8*delay(x,n,1)+0.1*delay(x,n,1)*delay(x,n,1);
end
figure
plot(1:3*N,y)
xlabel('n')
ylabel('y(n)')
title(fun)
%%%%%%%%%%%%%Configure noise%%%%%%%%%%%%%%%%
    P = 0
    if P == 0
        %%noise free
        e = zeros(1,3*N);
    elseif P == -1
        %%uniform noise
        rng(0);
        e = rand(1,3*N);
    else
        %%Gaussian noise
        sigma_e = sqrt(P/100*var(y)); 
        e = normrnd(0,sigma_e,[1,3*N]);
    end
    y_e = y+e;
    idealMSE = var(e)/var(y_e)*100;
    [a_opt pdelays_opt] = main(x,y,y_e,P,N,fun,idealMSE);
    disp(a_opt)
    disp(pdelays_opt)

function [a_opt pdelays_opt] = main(x,y,y_e,P,N,fun,idealMSE)
%% test for various delays and orders
tic;
min_mse = 99999;
for L = 1:6
   for K = 1:6
       N_0 = max(L,K);
       for order = 1:3
           [candidates delays] = genP(x(1,1:N),y_e(1,1:N),L,K,order,N);
           [a pdelays] = FOS(y_e(1,1:N),candidates,delays,N,N_0);
           [temp mse] = calcMSE(a,pdelays,x(N+1:2*N),y_e(N+1:2*N),N);
           if mse<min_mse
              min_mse = mse;
              a_opt = a;
              pdelays_opt = pdelays;
              K_opt = K;
              L_opt = L;
              order_opt = order;
           end
       end
   end
end
exetime = toc;
[y_1 mset] = calcMSE(a_opt,pdelays_opt,x(1:3*N),y(1:3*N),3*N);

disp(mset)
figure
plot(1:3*N,y)
hold on

plot(1:3*N,y_1)
xlabel('n')
ylabel('y(n)')
legend('y(n)','y_1(n)')
if P == 0
    title([fun,' without noise,',newline,'MSE=',num2str(mset),'%,runtime=',num2str(exetime),'s'])
elseif P == -1
    title([fun,' with uniform noise,',newline,'MSE=',num2str(mset),'%,runtime=',num2str(exetime),'s'])
else
    title([fun,' with Gaussian noise,',newline,'P=',num2str(P),'%,MSE=',num2str(mset),'%,runtime=',num2str(exetime),'s'])
end
hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [candidates tdelay] = genP(x,y,L,K,torder,N) %L= x delay, K= y delay
candidates = [];
tdelay = [];
for order=1:torder
    for o=0:order
        n_x = o;
        n_y = order-n_x;
        xdelay = combsrep(0:L,n_x);
        ydelay = combsrep(1:K,n_y);
        delays = cb(xdelay,ydelay,torder);
        dsize = size(delays);
        for i = 1:dsize(1)
            p = ones(1,N);
            for j = 1:dsize(2)-2
                if(isnan(delays(i,j)))
                    break;
                end
                if j<=n_x
                    for n = 1:N
                        p(n)=p(n)*delay(x,n,delays(i,j));
                    end
                else
                    for n = 1:N
                        p(n)=p(n)*delay(y,n,delays(i,j));
                    end
                end
            end
            candidates = [candidates;p];
        end
        tdelay = [tdelay;delays];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a pdelays] = FOS(y,candidates,delays,N,N_0)
%% initialization
csize = size(candidates);
p = ones(1,csize(2));              %p_0(n) = 1
dsize = size(delays);
pdelays = zeros(1,dsize(2));         %all 0 to represent constant term
Ncan = csize(1);                    %number of candidates remaining
trueg = mean(y);                    %g_0 = mean(y(n)w_0(n))/mean(w_0^2(n)) = mean(y(n))
M = 2;
Qsum = 0;
while true
    Q = 0;
    Qmax = 0;
    select = 0;
    g = trueg;
    for index = 1:Ncan
        D = 1;                     %D(0,0) = 1
        alpha = 0;
        p(M,:) = candidates(index,:);
        %% calc D and alpha
        for m=1+1:M
            D(m,0+1) = mean(p(m,:)); %D(m,0)=time average of p_m(n)%
            for r=(0+1):(m-1)
                alpha(m,r) = D(m,r)/D(r,r);
                i=(0+1):(r-1+1);
                D(m,r+1) = mean(p(m,:).*p(r+1,:)) - sum(alpha(r+1,i).*D(m,i));
            end
        end
        
        %% calc C          C_m= time average of y(n)w_m(n)
        C = mean(y);   %C_0=time average of y(n)w_0(n), where w_0(n)=p_0(n)=1
        for m=(1+1):M
            r = (0+1):(m-1);
            C(m) = mean(y.*p(m,:))-sum(alpha(m,r).*C(r));
        end
        
        %% find term with max Q and calc corresponding a
        g(M) = C(M)/D(M,M);
        Q = (g(M)^2)*D(M,M);
        if Q > Qmax
            Qmax = Q;
            select = index;
            truealpha = alpha;
            trueg = g;
        end
    end
    
    
    if (Qmax < 4*(mean(y(N_0:N).^2)-Qsum)/(N-N_0+1)) || select == 0
        %% terminate FOS
            M=M-1;
        for m = 0+1:M
            %%update corresponding v
            v(m) = 1;                          %v_m = 1;
            for i = m+1:M
                r = m:i-1;
                v(i) = -sum(truealpha(i,r).*v(r));
            end
            %%update corresponding a
            a(m) = sum(trueg(m:M).*v(m:M));
        end
        break;
    else
        %% FOS continues, need to update selected terms and remaining candidates
        p(M,:) = candidates(select,:);
        pdelays = [pdelays;delays(select,:)];   
        candidates(select,:) = [];
        delays(select,:) = [];
        M = M+1;
        Qsum = Qsum + Qmax;
        % update number of candidates remaining
        csize = size(candidates);
        Ncan = csize(1);
    end
end
end

%%%%%%%%%%%%%%%%%Util%%%%%%%%%%%%%%%%%%%%%%
function [y_1 mse] = calcMSE(a,pdelays,x,y,N)
temp = size(a);
NumTerm = temp(2);
y_1 = zeros(1,N);
for n = 1:N
    for i = 1:NumTerm
        y_1(n)=y_1(n)+a(i)*evalterm(x,y_1,n,pdelays(i,:));
    end
end
mse = mean(round(y_1-y,12).^2)/var(y)*100;   %I round it to 12 decimal since if I dont do this,for some reason even with exact same elements, y_1 ~= y. This not only fix this issue, but also improve my model selections
end

function val = delay(var,n,lag)
if n <= lag
    val = 0;
else
    val = var(n-lag);
end
end

function combs = cb(x,y,length)
xsize = size(x);
ysize = size(y);
ilength = length - xsize(2) - ysize(2);
if isempty(x)
    fill = NaN(ysize(1),ilength);
    for i=1:ysize(1)
        indicator(i,1)=0;
        indicator(i,2)=ysize(2);
    end
    combs = [y,fill,indicator];
elseif isempty(y)
    fill = NaN(xsize(1),ilength);
    for i=1:xsize(1)
        indicator(i,2)=0;
        indicator(i,1)=xsize(2);
    end
    combs = [x,fill,indicator];
else
    temp = 1;
    for ix = 1:xsize(1)
        for iy = 1:ysize(1)
            for j = 1:xsize(2)
                combs(temp,j) = x(ix,j);
            end
            for j = 1:ysize(2)
                combs(temp,j+xsize(2)) = y(iy,j);
            end
            for j = 1:ilength
                combs(temp,xsize(2)+ysize(2)+j) = NaN;
            end
            combs(temp,xsize(2)+ysize(2)+ilength+1) = xsize(2);
            combs(temp,xsize(2)+ysize(2)+ilength+2) = ysize(2);
            temp = temp+1;
        end
    end
end
end

function prod = evalterm(x,y,n,pdelays)
dsize = size(pdelays);
dcoln = dsize(2);
prod = 1;
numX = pdelays(1,dcoln-1);
numY = pdelays(1,dcoln);
if (numX == 0) && (numY == 0)
    prod = 1;
else
    if numX ~= 0
        for j = 1:numX
            prod = prod * delay(x,n,pdelays(1,j));
        end
    end
    if numY ~= 0
        for j = numX+1:numX+numY
            prod = prod * delay(y,n,pdelays(1,j));
        end
    end
end

end
