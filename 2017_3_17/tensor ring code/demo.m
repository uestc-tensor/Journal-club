%% experiment 1 function 1

clear all;

% generate fuction y
x=-1:1.45e-06:1;
y=(x+1).*(sin(100.*(x+1).^2));
plot(y);
y=y(1,1:1048576);
Y=reshape(y,4,4,4,4,4,4,4,4,4,4);%construct 10-dimensinal tensor
dim=size(Y);
d=length(dim);
tot=3e-4;
%% approximation using trsvd
[Zsvd,r]=trsvd(Y,tot);
average=sum(r)/length(r);

%% using TR-alsmaxiter=100;
tot=3e-4;
r(10+1)=r(1);
[Zals]=trals(Y,r,maxiter,tot);
[Y_pre]=constract_X(Zals,r,d,dim); %compute approximation of Y
Y_pre=reshape(Y_pre,1,[]);
 plot(Y_pre);

%% trbls
tot=9e-4;
[Zbals,r]=trbals(Y,maxiter,tot);
[Y_pre]=constract_X(Zbals,r,d,dim); %compute approximation of Y
Y_pre=reshape(Y_pre,1,[]);
 plot(Y_pre);

