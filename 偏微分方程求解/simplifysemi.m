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
P1t=ones(128,128,36)
P2t=ones(128,128,36);
P3t=ones(128,128,36);
P1tF=ones(128,128,36);
P2tF=ones(128,128,36);
P3tF=ones(128,128,36);
% x1=ones(123,128,36);
% x2=ones(128,128,36);
% x3=ones(128,128,36);
% for k=1:128
%     x1(k,:,:)=k;
% end
% for k=1:128
%     x2(:,k,:)=k;
% end
% for k=1:36
%     x3(:,:,k)=(k-4);
% end
xx=1:128;
yy=1:128;
zz=-3:32;
[x1,x2,x3]=meshgrid(xx,yy,zz);
bian1=x1;
bian2=x2;
bian3=x3;
P1t=P1;
P2t=P2;
P3t=P3;
P1fourier=fftn(P1);
P2fourier=fftn(P2);
P3fourier=fftn(P3);
P1tF=P1fourier;
P2tF=P2fourier;
P3tF=P3fourier;

fo=fsuo(P1t,P2t,P3t);
foF=fftn(fo);
bian=sqrt(bian1.^2+bian2.^2+bian3.^2);

P1tF=(P1tF+detat*foF)./(detat*bian.^2+1);
P2tF=(P2tF+detat*foF)./(detat*bian.^2+1);
P3tF=(P3tF+detat*foF)./(detat*bian.^2+1);
F1t=P1t;
F2t=P2t;
F3t=P3t;
P1t=real(ifftn(P1tF));
P2t=real(ifftn(P2tF));
P3t=real(ifftn(P3tF));
k=k+1;
while k<50
   fo=fsuo(P1t,P2t,P3t);
   foo=fsuo(P1t,P2t,P3t);
foF=fftn(fo);
fooF=fftn(foo);
P1tF=(4*P1tF-P1tF+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
P2tF=(4*P2tF-P2tF+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
P3tF=(4*P3tF-P3tF+2*detat*(2*foF-fooF))./(3+2*detat*bian.^2);
F1t=P1t;
F2t=P2t;
F3t=P3t;
P1t=real(ifftn(P1tF));
P2t=real(ifftn(P2tF));
P3t=real(ifftn(P3tF));
P3t(1,1,1)
    k=k+1;
end
[x,y,z]=meshgrid(1:1:128,1:1:128,1:1:36);
 PTR1=((P1t)/(max(max(max(P1t)))-min(min(min(P1t)))))*225;
 figure (1)
 slice(x,y,z,PTR1,[1,128],[1,128],[1,36]);
PTR2=((P2t)/(max(max(max(P2t)))-min(min(min(P2t)))))*225;
figure (2)
slice(x,y,z,PTR2,[1,128],[1,128],[1,36]);
PTR3=((P3t)/(max(max(max(P3t)))-min(min(min(P3t)))))*225;
figure (3)
slice(x,y,z,PTR3,[1,128],[1,128],[1,36]);