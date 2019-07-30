detat=0.000001;
yita=ones(1,128);
for i=1:128
    if(i<65)
        yita(i)=1;
    else
        yita(i)=-1;
    end
end
x1=1:128;
bian1=x1;
P1=yita;
P1t(1,:)=P1;
P1fourier=fftn(P1);
P1tF(1,:)=P1fourier;
k=1;
fo=ted(P1);
foF=fftn(fo);
bian=sqrt(bian1.^2);
P1tF(k+1,:)=(P1tF(k,:)+detat*foF)./(detat*bian.^2+1);
P1t(k+1,:)=real(ifftn(P1tF(k+1,:)));
k=k+1;
c(k)=abs(mean(abs(P1t(k,:)))-mean(abs(P1t(k-1,:))));
while c(k)>mean(abs(P1t(k,:)))/1000000
    fo=ted(P1t(k,:));
    foo=ted(P1t(k-1,:));
    foF=fftn(fo);
    fooF=fftn(foo);
    P1tF(k+1,:)=(4*P1tF(k,:)-P1tF(k-1,:)+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
    P1t(k+1,:)=real(ifftn(P1tF(k+1,:)));
    P1t(k+1,33)
   k=k+1; 
   c(k)= abs(mean(abs(P1t(k,:)))-mean(abs(P1t(k-1,:))));
end
plot(x1,P1t(k-1,:));