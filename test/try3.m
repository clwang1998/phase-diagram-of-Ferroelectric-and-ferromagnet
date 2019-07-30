detat=0.0001;
yita=ones(1,128);
n=128;
yita(1:n/2)=1;
yita(n/2+1:128)=-1;
x1=1:128;
P1=yita;
P1t(1,:)=P1;
P1fourier=fft(P1);
P1tF(1,:)=P1fourier;
k=0;
while k<10000
    k=k+1;
    fo=ted(P1t(k,:));
    foF=fft(fo);
    for l=1:128
    P1tF(k+1,l)=P1tF(k,l)+(foF(l)-P1tF(k,l)*abs(P1tF(k,l))^2)*detat;
    P1t(k+1,:)=(ifft(P1tF(k+1,:)));
    end
    P1t(k,1)
   
end