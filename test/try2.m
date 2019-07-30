clc;clear;
P1=randn(128,128,20);
P2=randn(128,128,20);
P3=randn(128,128,20);
dx=1e-9;

p11=diff(P1,1,1)/dx;P11=cat(1,p11(1,:,:),p11);
p22=diff(P2,1,2)/dx;P22=cat(2,p22(:,1,:),p22);
p33=diff(P3,1,3)/dx;P33=cat(3,p33(:,:,1),p33);

yibu0=8.85*1e-12;
k11=190;k22=190;k33=190;
s=(P11+P22+P33)/yibu0;

fi=zeros(128,128,20);
fiold=fi+1;
finew=fi;
k=1;
while k<10000
    if max(max(max(abs(fi-fiold))))>=1e-4
        finew(2:127,2:127,2:19)=(k11*(fi(3:128,2:127,2:19)+fi(1:126,2:127,2:19))+...
            k22*(fi(2:127,3:128,2:19)+fi(2:127,1:126,2:19))+...
            k33*(fi(2:127,2:127,3:20)+fi(2:127,2:127,1:18))-dx^2*s(2:127,2:127,2:19))/2/(k11+k22+k33);
        fiold=fi;
        fi=finew;
    else
        k=10000;
    end
    k=k+1;
end

[E1,E2,E3]=gradient(fi);
E1=-E1;E2=-E2;E3=-E3;
felec=-0.5*(E1.*P1(:,:,1:20)+E2.*P2(:,:,1:20)+E3.*P3(:,:,1:20));%%%%
felec=cat(3,felec,zeros(128,128,16));
Felec=sum(sum(sum(felec*(dx)^3)));
 