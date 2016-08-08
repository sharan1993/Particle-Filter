clear;
data=load('magnets-data.txt');

%given values
sig_a=0.0625; 
xm1=-10;
xm2=10;
sig_n=0.003906;
sig_m=4.0;

%Initialization
T=1;                           %sampling rate
N=1000;                          %No. of particles
Yt=data(:,3);
for i=1:N
    x(1,i)=0;                   %Position Variable
    x_dot(1,i)=0;               %Velocity variable
    y(1,i)=0;                %output..sensor reading
    w(1,i)=1/N;                 %initial weights
    pdf(1,i)=0;                      %probability distribution function
end

%state transition equation
for t=2:length(data(:,1))
    Exp=0;CV=0;
for i=1:N
   x(t,i)=x(t-1,i)+x_dot(t-1,i)*T;
   a_t=normrnd(0,(sig_a));
   if (x(t-1,i)<-20)
        x_dot(t,i)=2;
   elseif (-20<=x(t-1,i) && x(t-1,i)<0)
        x_dot(t,i)=x_dot(t-1,i)+abs(a_t);
   elseif (0<=x(t-1,i) && x(t-1,i)<=20)
        x_dot(t,i)=x_dot(t-1,i)-abs(a_t);
   elseif (x(t-1,i)>20)
        x_dot(t,i)=-2;
   end
   n_t=normrnd(0,sig_n^2);
    %y(t,i)=normpdf(x(t,i),xm1,sig_m)+normpdf(x(t,i),xm2,sig_m);
    y(t,i)=normpdf(x(t,i),xm1,sig_m)+normpdf(x(t,i),xm2,sig_m);
    %im(t,i)=normpdf(x(t,i),xm1,sig_m)+normpdf(x(t,i),xm2,sig_m);
    pdf(t,i)=normpdf(y(t,i),data(t,3),sig_n);
    %pdf(t-1,i)=normpdf(im(t-1,i),y(t-1,i),sig_n);
    w(t,i)=w(t-1,i)*pdf(t,i);
end

for i=1:N
    wn(t,i)=w(t,i)/sum(w(t,:));
    Exp=Exp+x(t-1,i)*wn(t,i);
end
w=wn;
    E(t)=Exp;
    CV=0;
    for i=1:N
        CV=CV+(N*w(t,i)-1)^2;
    end
    CV=CV/N;
    ESS=N/(1+CV);
if(ESS<0.5*N) 
    asd=[];
    Q=cumsum(w(t,:));
    t1=rand(N+1);
    t2=sort(t1);
    t2(N+1)=1.0;
    j=1;k=1;
    while (j<=N)
        if (t2(j)<Q(k))
            asd(j)=k;
            j=j+1;
        else 
            k=k+1;
        end
    end
for i=1:N
    x(t,i)=x(t,asd(i));
    x_dot(t,i)=x_dot(t,asd(i));
    w(t,i)=1/N;

end
end
end
