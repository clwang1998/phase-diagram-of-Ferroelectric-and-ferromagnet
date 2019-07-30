detat=0.00001;
P3=ones(128,128,36)*1e-6;
kk=1;
while kk<64*64*10
    P3(randperm(128,1),randperm(128,1),randperm(20,1))=P3(randperm(128,1),randperm(128,1),randperm(20,1))+randn*1e-6;
    kk=kk+1;
end
P3(:,:,21:36)=0;
P1=zeros(128,128,36);
P2=zeros(128,128,36);
P1t=ones(128,128,36,100);
P2t=ones(128,128,36,100);
P3t=ones(128,128,36,100);
P1tF=ones(128,128,36,100);
P2tF=ones(128,128,36,100);
P3tF=ones(128,128,36,100);
x1=ones(123,128,36);
x2=ones(128,128,36);
x3=ones(128,128,36);
for k=1:128
    x1(k,:,:)=k;
end
for k=1:128
    x2(:,k,:)=k;
end
for k=1:36
    x3(:,:,k)=(k-4);
end
bian1=x1;
bian2=x2;
bian3=x3;
P1t(:,:,:,1)=P1;
P2t(:,:,:,1)=P2;
P3t(:,:,:,1)=P3;
P1fourier=fftn(P1);
P2fourier=fftn(P2);
P3fourier=fftn(P3);
P1tF(:,:,:,1)=P1fourier;
P2tF(:,:,:,1)=P2fourier;
P3tF(:,:,:,1)=P3fourier;
k=1;
fo=fsuo(P1t(:,:,:,k),P2t(:,:,:,k),P3t(:,:,:,k));
foF=fftn(fo);
bian=sqrt(bian1.^2+bian2.^2+bian3.^2);

P1tF(:,:,:,k+1)=(P1tF(:,:,:,k)+detat*foF)./(detat*bian.^2+1);
P2tF(:,:,:,k+1)=(P2tF(:,:,:,k)+detat*foF)./(detat*bian.^2+1);
P3tF(:,:,:,k+1)=(P3tF(:,:,:,k)+detat*foF)./(detat*bian.^2+1);
P1t(:,:,:,k+1)=real(ifftn(P1tF(:,:,:,k+1)));
P2t(:,:,:,k+1)=real(ifftn(P2tF(:,:,:,k+1)));
P3t(:,:,:,k+1)=real(ifftn(P3tF(:,:,:,k+1)));
k=k+1;
while k<100
   fo=fsuo(P1t(:,:,:,k),P2t(:,:,:,k),P3t(:,:,:,k));
   foo=fsuo(P1t(:,:,:,k-1),P2t(:,:,:,k-1),P3t(:,:,:,k-1));
foF=fftn(fo);
fooF=fftn(foo);
P1tF(:,:,:,k+1)=(4*P1tF(:,:,:,k)-P1tF(:,:,:,k-1)+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
P2tF(:,:,:,k+1)=(4*P2tF(:,:,:,k)-P2tF(:,:,:,k-1)+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
P3tF(:,:,:,k+1)=(4*P3tF(:,:,:,k)-P3tF(:,:,:,k-1)+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
P1t(:,:,:,k+1)=real(ifftn(P1tF(:,:,:,k+1)));
P2t(:,:,:,k+1)=real(ifftn(P2tF(:,:,:,k+1)));
P3t(:,:,:,k+1)=real(ifftn(P3tF(:,:,:,k+1)));
P3t(1,1,1,k+1)
    k=k+1;
end
[x,y,z]=meshgrid(1:1:128,1:1:128,1:1:36);
 PTR1=((P1t(:,:,:,k))/(max(max(max(P1t(:,:,:,k))))-min(min(min(P1t(:,:,:,k))))))*225;
 figure (1)
 slice(x,y,z,PTR1,[1,128],[1,128],[1,36]);
PTR2=((P2t(:,:,:,k))/(max(max(max(P2t(:,:,:,k))))-min(min(min(P2t(:,:,:,k))))))*225;
figure (2)
slice(x,y,z,PTR2,[1,128],[1,128],[1,36]);
PTR3=((P3t(:,:,:,k))/(max(max(max(P3t(:,:,:,k))))-min(min(min(P3t(:,:,:,k))))))*225;
figure (3)
slice(x,y,z,PTR3,[1,128],[1,128],[1,36]);