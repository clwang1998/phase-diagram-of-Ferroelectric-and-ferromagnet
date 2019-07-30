function f = f(P1,P2,P3 )
dx0=1e-9;
e0=0;  T=25;            %%%%
dx=(1+e0)*dx0;
G11=dx^2*abs(4.124*(25-115)*1e5);


%F_bulk
a1=4.124*(T-115)*1e5;a11=-1.950*1e9;a12=7.974*1e8;a111=1.294*1e9;a112=-1.950*1e9;a123=-2.5*1e9;
a1111=3.860*1e10;a1112=2.529*1e10;a1122=1.637*1e10;a1123=1.367*1e10;
fbulk=a1*(P1.^2+P2.^2+P3.^2)+a11*(P1.^4+P2.^4+P3.^4)+a12*(P1.^2.*P2.^2+P3.^2.*P2.^2+P1.^2.*P3.^2)+a111*(P1.^6+P2.^6+P3.^6)+...
    a112*(P1.^2.*(P2.^4+P3.^4)+P2.^2.*(P1.^4+P3.^4)+P3.^2.*(P2.^4+P1.^4))+a123*P1.^2.*P2.^2.*P3.^2+...
    a1111*(P1.^8+P2.^8+P3.^8)+a1112*(P1.^6.*(P2.^2+P3.^2)+P2.^6.*(P1.^2+P3.^2)+P3.^6.*(P2.^2+P1.^2))+...
    a1122*(P1.^4.*P2.^4+P3.^4.*P2.^4+P1.^4.*P3.^4)+a1123*(P1.^4.*P2.^2.*P3.^2+P1.^2.*P2.^4.*P3.^2+P1.^2.*P2.^2.*P3.^4);
Fbulk=sum(sum(sum(fbulk*(dx)^3)));

%F_wall
p11=diff(P1,1,1)/dx;P11=cat(1,p11(1,:,:),p11);
p12=diff(P1,1,2)/dx;P12=cat(2,p12(:,1,:),p12);
p13=diff(P1,1,3)/dx;P13=cat(3,p13(:,:,1),p13);
p21=diff(P2,1,1)/dx;P21=cat(1,p21(1,:,:),p21);
p22=diff(P2,1,2)/dx;P22=cat(2,p22(:,1,:),p22);
p23=diff(P2,1,3)/dx;P23=cat(3,p23(:,:,1),p23);
p31=diff(P3,1,1)/dx;P31=cat(1,p31(1,:,:),p31);
p32=diff(P3,1,2)/dx;P32=cat(2,p32(:,1,:),p32);
p33=diff(P3,1,3)/dx;P33=cat(3,p33(:,:,1),p33);
fwall=0.5*G11*(P11.^2+P12.^2+P13.^2+P21.^2+P22.^2+P23.^2+P31.^2+P32.^2+P33.^2);
Fwall=sum(sum(sum(fwall*(dx)^3)));

%F_elec
P1bian=1/(2*pi)^3*fftn(P1(:,:,1:20));
P2bian=1/(2*pi)^3*fftn(P2(:,:,1:20));
P3bian=1/(2*pi)^3*fftn(P3(:,:,1:20));

x1=ones(128,128,20)*dx;
x2=ones(128,128,20)*dx;
x3=ones(128,128,20)*dx;
for i=1:128
    x1(i,:,:)=i*dx;
end
for i=1:128
    x2(:,i,:)=i*dx;
end
for i=1:20
    x3(:,:,i)=(21-i)*dx;
end
bian1=x1/dx;
bian2=x2/dx;
bian3=ones(128,128,20);
for i=1:20
    bian3(:,:,i)=i;
end

yibu0=8.85*1e-12;
k11=190;k22=190;k33=190;
fiAbian=-i*(bian1.*P1bian+bian2.*P2bian+bian3.*P3bian)./(yibu0*(k11*bian1.^2+k22*bian2.^2+k33*bian3.^2));
fiA=ifftn(fiAbian);

fi1=0;fi2=0;   %短路电学边界条件
D1=fi1-fiA(:,:,20);
D2=fi2-fiA(:,:,1);
W1=exp(20*dx*((k11*bian1(:,:,1).^2+k22*bian2(:,:,1).^2).^0.5));
W2=1./W1;
C1=(D1.*W1-D2)./(W1-W2);
C2=D1-C1;

C1a=C1;
for i=1:19
    C1=cat(3,C1,C1a);
end
C2a=C2;
for i=1:19
    C2=cat(3,C2,C2a);
end
fiBbian=C1.*exp(x3.*(((k11*bian1.^2+k22*bian2.^2)/k33).^0.5))+...
    C2.*exp(-x3.*(((k11*bian1.^ 2+k22*bian2.^2)/k33).^0.5));
fiB=ifft(ifft(fiBbian,[],1),[],2);

fi=real(fiA+fiB);
[E1,E2,E3]=gradient(fi);
E1=-E1;E2=-E2;E3=-E3;
felec=-0.5*(E1.*P1(:,:,1:20)+E2.*P2(:,:,1:20)+E3.*P3(:,:,1:20));%%%%
felec=cat(3,felec,zeros(128,128,16));
Felec=sum(sum(sum(felec*(dx)^3)));

% F_elas
C11=1.78*10^11;
C12=0.964*10^11;
C44=1.22*10^11;
Q11=0.1;
Q12=-0.034;
Q44=0.029;


%应力场
yibu011=Q11*P1.^2+Q12*(P2.^2+P3.^2);
yibu022=Q11*P2.^2+Q12*(P1.^2+P3.^2);
yibu033=Q11*P3.^2+Q12*(P2.^2+P1.^2);
yibu023=2*Q44*P2.*P3;
yibu013=2*Q44*P1.*P3;
yibu012=2*Q44*P1.*P2;
yibu011B=1/(2*pi)^3*fftn(yibu011(:,:,1:32));
yibu022B=1/(2*pi)^3*fftn(yibu022(:,:,1:32));
yibu033B=1/(2*pi)^3*fftn(yibu033(:,:,1:32));
yibu023B=1/(2*pi)^3*fftn(yibu023(:,:,1:32));
yibu013B=1/(2*pi)^3*fftn(yibu013(:,:,1:32));
yibu012B=1/(2*pi)^3*fftn(yibu012(:,:,1:32));
yibu032B=yibu023B;
yibu031B=yibu013B;
yibu021B=yibu012B;

x1=ones(128,128,32)*dx;
x2=ones(128,128,32)*dx;
x3=ones(128,128,32)*dx;
for i=1:128
    x1(i,:,:)=i*dx;
end
for i=1:128
    x2(:,i,:)=i*dx;
end
for i=1:32
    x3(:,:,i)=(21-i)*dx;
end
bian1=x1/dx;
bian2=x2/dx;
bian3=ones(128,128,32);
for i=1:32
    bian3(:,:,i)=i;
end



hhhj1=(C11*yibu011B+C12*yibu022B+C12*yibu033B).*bian1+...
    (C44*yibu012B+C44*yibu021B).*bian2+(C44*yibu013B+C44*yibu031B).*bian3;
hhhj2=(C44*yibu012B+C44*yibu021B).*bian1+(C12*yibu011B+C11*yibu022B+C12*yibu033B).*bian2+...
    2*C44*yibu023B.*bian3;
hhhj3=2*C44*yibu013B.*bian1+2*C44*yibu023B.*bian2+C12*yibu011B.*bian3+C12*yibu022B.*bian3+...
    C11*yibu033B.*bian3;
for i=1:128
    for j=1:128
        for k=1:32
            G(1,1)=C11*bian1(i,j,k)^2+C44*bian2(i,j,k)^2+C44*bian3(i,j,k)^2;
            G(1,2)=(C12+C44)*bian1(i,j,k)*bian2(i,j,k);
            G(1,3)=(C12+C44)*bian1(i,j,k)*bian3(i,j,k);
            G(2,1)=G(1,2);
            G(2,2)=C11*bian2(i,j,k)^2+C44*bian1(i,j,k)^2+C44*bian3(i,j,k)^2;
            G(2,3)=(C12+C44)*bian2(i,j,k)*bian3(i,j,k);
            G(3,1)=G(1,3);
            G(3,2)=G(2,3);
            G(3,3)=C11*bian3(i,j,k)^2+C44*bian1(i,j,k)^2+C44*bian2(i,j,k)^2;
            g=inv(G);
            uA1bian(i,j,k)=-i*(g(1,1)*hhhj1(i,j,k)+g(1,2)*hhhj2(i,j,k)+g(1,3)*hhhj3(i,j,k));
            uA2bian(i,j,k)=-i*(g(2,1)*hhhj1(i,j,k)+g(2,2)*hhhj2(i,j,k)+g(2,3)*hhhj3(i,j,k));
            uA3bian(i,j,k)=-i*(g(3,1)*hhhj1(i,j,k)+g(3,2)*hhhj2(i,j,k)+g(3,3)*hhhj3(i,j,k));
        end
    end
end
uA1=ifftn(uA1bian);
uA2=ifftn(uA2bian);
uA3=ifftn(uA3bian);

%uB
u=100;
v=0.1;
bian1=bian1(:,:,1);
bian2=bian2(:,:,1);
bian=sqrt(bian1.^2+bian2.^2);

uA1bian=1/(2*pi)^2*fft(fft(uA1,[],1),[],2);
uA2bian=1/(2*pi)^2*fft(fft(uA2,[],1),[],2);
uA3bian=1/(2*pi)^2*fft(fft(uA3,[],1),[],2);
uA33bian=diff(uA1bian,1,1);uA33hf=uA33bian(:,:,1);
uA13bian=diff(uA1bian,1,3);uA13hf=uA13bian(:,:,1);
uA23bian=diff(uA2bian,1,3);uA23hf=uA23bian(:,:,1);
yibu011B=1/(2*pi)^2*fft(fft(yibu011,[],1),[],2);
yibu012B=1/(2*pi)^2*fft(fft(yibu012,[],1),[],2);
yibu013B=1/(2*pi)^2*fft(fft(yibu013,[],1),[],2);
yibu022B=1/(2*pi)^2*fft(fft(yibu022,[],1),[],2);
yibu033B=1/(2*pi)^2*fft(fft(yibu033,[],1),[],2);
yibu023B=1/(2*pi)^2*fft(fft(yibu023,[],1),[],2);
D1=-(i*bian1.*uA3bian(:,:,1)+uA13hf-2*yibu013B(:,:,1));
D2=-(i*bian2.*uA3bian(:,:,1)+uA23hf-2*yibu023B(:,:,1));
D3=-(C12*i*bian1.*uA1bian(:,:,1)+C12*i*bian2.*uA2bian(:,:,1)+C11*uA3bian(:,:,1))+...
    (C12*yibu011B(:,:,1)+C12*yibu022B(:,:,1)+C11*yibu033B(:,:,1));
D4=-uA1bian(:,:,32);
D5=-uA2bian(:,:,32);
D6=-uA3bian(:,:,32);

a11=-i/u*bian2./bian1;a12=i/u*ones(128,128);a13=zeros(128,128);
a21=-i/(4*u*v)*bian1./bian;a22=i/(4*u*v)*bian2./bian;a23=-1/(4*u*v)*ones(128,128);
a31=(1-2*v)/(2*u*v)*bian./bian1;a32=zeros(128,128);a33=i/(4*u*v)*ones(128,128);
a41=conj(a11);a42=conj(a12);a43=conj(a13);
a51=conj(a21);a52=conj(a22);a53=conj(a23);
a61=conj(a31);a62=conj(a32);a63=conj(a33);

hf=20*dx;
cu11f=a11.*exp(-bian*hf);
cu12f=a12.*exp(-bian*hf);
cu13f=a13.*exp(-bian*hf);
cu21f=a21.*exp(-bian*hf);
cu22f=a22.*exp(-bian*hf);
cu23f=a23.*exp(-bian*hf);
cu31f=a31.*exp(-bian*hf)+1i*hf*a21.*exp(-bian*hf);
cu32f=a32.*exp(-bian*hf)+1i*hf*a22.*exp(-bian*hf);
cu33f=a33.*exp(-bian*hf)+1i*hf*a23.*exp(-bian*hf);
cu41f=a41.*exp(bian*hf);
cu42f=a42.*exp(bian*hf);
cu43f=a43.*exp(bian*hf);
cu51f=a51.*exp(bian*hf);
cu52f=a52.*exp(bian*hf);
cu53f=a53.*exp(bian*hf);
cu61f=a61.*exp(bian*hf)+1i*hf*a51.*exp(bian*hf);
cu62f=a62.*exp(bian*hf)+1i*hf*a52.*exp(bian*hf);
cu63f=a63.*exp(bian*hf)+1i*hf*a53.*exp(bian*hf);
eu11=-bian.*a11.*exp(-bian*hf);
eu12=-bian.*a12.*exp(-bian*hf);
eu13=-bian.*a13.*exp(-bian*hf);
eu21=-bian.*a21.*exp(-bian*hf);
eu22=-bian.*a22.*exp(-bian*hf);
eu23=-bian.*a23.*exp(-bian*hf);
eu31=-bian.*a31.*exp(-bian*hf)+1i*a21.*exp(-bian*hf)-1i*hf*bian.*a21.*exp(-bian*hf);
eu32=-bian.*a32.*exp(-bian*hf)+1i*a22.*exp(-bian*hf)-1i*hf*bian.*a22.*exp(-bian*hf);
eu33=-bian.*a33.*exp(-bian*hf)+1i*a23.*exp(-bian*hf)-1i*hf*bian.*a23.*exp(-bian*hf);
eu41=bian.*a41.*exp(bian*hf);
eu42=bian.*a42.*exp(bian*hf);
eu43=bian.*a43.*exp(bian*hf);
eu51=bian.*a51.*exp(bian*hf);
eu52=bian.*a52.*exp(bian*hf);
eu53=bian.*a53.*exp(bian*hf);
eu61=bian.*a51.*exp(bian*hf)+1i*a51.*exp(bian*hf)+1i*hf*bian.*a51.*exp(bian*hf);
eu62=bian.*a52.*exp(bian*hf)+1i*a52.*exp(bian*hf)+1i*hf*bian.*a52.*exp(bian*hf);
eu63=bian.*a53.*exp(bian*hf)+1i*a53.*exp(bian*hf)+1i*hf*bian.*a53.*exp(bian*hf);

hs=12*dx;
hss=-hs;
cu11s=a11.*exp(-bian*hss);
cu12s=a12.*exp(-bian*hss);
cu13s=a13.*exp(-bian*hss);
cu21s=a21.*exp(-bian*hss);
cu22s=a22.*exp(-bian*hss);
cu23s=a23.*exp(-bian*hss);
cu31s=a31.*exp(-bian*hss)+1i*hss*a21.*exp(-bian*hss);
cu32s=a32.*exp(-bian*hss)+1i*hss*a22.*exp(-bian*hss);
cu33s=a33.*exp(-bian*hss)+1i*hss*a23.*exp(-bian*hss);
cu41s=a41.*exp(bian*hss);
cu42s=a42.*exp(bian*hss);
cu43s=a43.*exp(bian*hss);
cu51s=a51.*exp(bian*hss);
cu52s=a52.*exp(bian*hss);
cu53s=a53.*exp(bian*hss);
cu61s=a61.*exp(bian*hss)+1i*hss*a51.*exp(bian*hss);
cu62s=a62.*exp(bian*hss)+1i*hss*a52.*exp(bian*hss);
cu63s=a63.*exp(bian*hss)+1i*hss*a53.*exp(bian*hss);

Q=ones(128,128,6);
for aa=1:128
    for bb=1:128
A=[1i*bian1(aa,bb)*cu13f(aa,bb)+eu11(aa,bb),1i*bian1(aa,bb)*cu23f(aa,bb)+eu21(aa,bb),1i*bian1(aa,bb)*cu33f(aa,bb)+eu31(aa,bb),1i*bian1(aa,bb)*cu43f(aa,bb)+eu41(aa,bb),1i*bian1(aa,bb)*cu53f(aa,bb)+eu51(aa,bb),1i*bian1(aa,bb)*cu63f(aa,bb)+eu61(aa,bb);
   1i*bian2(aa,bb)*cu13f(aa,bb)+eu12(aa,bb),1i*bian2(aa,bb)*cu23f(aa,bb)+eu22(aa,bb),1i*bian2(aa,bb)*cu33f(aa,bb)+eu32(aa,bb),1i*bian2(aa,bb)*cu43f(aa,bb)+eu42(aa,bb),1i*bian2(aa,bb)*cu53f(aa,bb)+eu52(aa,bb),1i*bian2(aa,bb)*cu63f(aa,bb)+eu62(aa,bb);
   1i*C12*bian1(aa,bb)*cu11f(aa,bb)+1i*C12*bian2(aa,bb)*cu12f(aa,bb)+C11*eu13(aa,bb),1i*C12*bian1(aa,bb)*cu21f(aa,bb)+1i*C12*bian2(aa,bb)*cu22f(aa,bb)+C11*eu23(aa,bb),1i*C12*bian1(aa,bb)*cu31f(aa,bb)+1i*C12*bian2(aa,bb)*cu32f(aa,bb)+C11*eu33(aa,bb),1i*C12*bian1(aa,bb)*cu41f(aa,bb)+1i*C12*bian2(aa,bb)*cu42f(aa,bb)+C11*eu43(aa,bb),1i*C12*bian1(aa,bb)*cu51f(aa,bb)+1i*C12*bian2(aa,bb)*cu52f(aa,bb)+C11*eu53(aa,bb),1i*C12*bian1(aa,bb)*cu61f(aa,bb)+1i*C12*bian2(aa,bb)*cu62f(aa,bb)+C11*eu63(aa,bb);
   cu11s(aa,bb),cu21s(aa,bb),cu31s(aa,bb),cu41s(aa,bb),cu51s(aa,bb),cu61s(aa,bb);
   cu12s(aa,bb),cu22s(aa,bb),cu32s(aa,bb),cu42s(aa,bb),cu52s(aa,bb),cu62s(aa,bb);
   cu13s(aa,bb),cu23s(aa,bb),cu33s(aa,bb),cu43s(aa,bb),cu53s(aa,bb),cu63s(aa,bb)];
B=[D1(aa,bb);D2(aa,bb);D3(aa,bb);D4(aa,bb);D5(aa,bb);D6(aa,bb)];
test=B\A;
Q(aa,bb,1)=test(1);
Q(aa,bb,2)=test(2);
Q(aa,bb,3)=test(3);
Q(aa,bb,4)=test(4);
Q(aa,bb,5)=test(5);
Q(aa,bb,6)=test(6);
    end
end
q1=Q(:,:,1);
q2=Q(:,:,2);
q3=Q(:,:,3);
q4=Q(:,:,4);
q5=Q(:,:,5);
q6=Q(:,:,6);

q1a=q1;
for i=1:31
    q1=cat(3,q1,q1a);
end
q2a=q2;
for i=1:31
    q2=cat(3,q2,q2a);
end
q3a=q3;
for i=1:31
    q3=cat(3,q3,q3a);
end
q4a=q4;
for i=1:31
    q4=cat(3,q4,q4a);
end
q5a=q5;
for i=1:31
    q5=cat(3,q5,q5a);
end
q6a=q6;
for i=1:31
    q6=cat(3,q6,q6a);
end

uB1bian=q1.*a11.*exp(-bian.*x3)+q2.*a21.*exp(-bian.*x3)+q3.*(a31.*exp(-bian.*x3)+i*bian.*x3.*a21.*exp(-bian.*x3))+...
    q4.*a41.*exp(bian.*x3)+q5.*a51.*exp(bian.*x3)+q6.*(a61.*exp(bian.*x3)+i*bian.*x3.*a51.*exp(bian.*x3));
uB2bian=q1.*a12.*exp(-bian.*x3)+q2.*a22.*exp(-bian.*x3)+q3.*(a32.*exp(-bian.*x3)+i*bian.*x3.*a22.*exp(-bian.*x3))...
       +q4.*a42.*exp(bian.*x3)+q5.*a52.*exp(bian.*x3)+q6.*(a62.*exp(bian.*x3)+i*bian.*x3.*a52.*exp(bian.*x3));
uB3bian=q1.*a13.*exp(-bian.*x3)+q2.*a23.*exp(-bian.*x3)+q3.*(a33.*exp(-bian.*x3)+i*bian.*x3.*a23.*exp(-bian.*x3))...
       +q4.*a43.*exp(bian.*x3)+q5.*a53.*exp(bian.*x3)+q6.*(a63.*exp(bian.*x3)+i*bian.*x3.*a53.*exp(bian.*x3));
uB1=ifft(ifft(uB1bian,[],1),[],2);
uB2=ifft(ifft(uB2bian,[],1),[],2);
uB3=ifft(ifft(uB3bian,[],1),[],2);

u1=real(uA1+uB1);
u2=real(uA2+uB2);
u3=real(uA3+uB3);
u11=diff(u1,1,1)/dx;u11=cat(1,u11(1,:,:),u11);
u12=diff(u1,1,2)/dx;u12=cat(2,u12(:,1,:),u12);
u13=diff(u1,1,3)/dx;u13=cat(3,u13(:,:,1),u13);
u21=diff(u2,1,1)/dx;u21=cat(1,u21(1,:,:),u21);
u22=diff(u2,1,2)/dx;u22=cat(2,u22(:,1,:),u22);
u23=diff(u2,1,3)/dx;u23=cat(3,u23(:,:,1),u23);
u31=diff(u3,1,1)/dx;u31=cat(1,u31(1,:,:),u31);
u32=diff(u3,1,2)/dx;u32=cat(2,u32(:,1,:),u32);
u33=diff(u3,1,3)/dx;u33=cat(3,u33(:,:,1),u33);
yin11=u11;
yin22=u22;
yin33=u33;
yin13=0.5*(u13+u31);
yin12=0.5*(u12+u21);
yin23=0.5*(u23+u32);
yin31=yin13;yin32=yin23;yin21=yin12;

yibu11=e0+yin11;
yibu22=e0+yin22;
yibu33=-C12*2/C11*e0+yin33;%%这里三个正应变有约束条件
yibu12=0+yin12;yibu21=0+yin21;
yibu13=0+yin13;yibu31=0+yin31;
yibu23=0+yin23;yibu32=0+yin32;


felas1=0.5*C11*(yibu11.^2+yibu22.^2+yibu33.^2)+C12*(yibu11.*yibu22+yibu22.*yibu33+yibu11.*yibu33)+...
    2*C44*(yibu13.^2+yibu12.^2+yibu23.^2);
Felas1=sum(sum(sum(felas1*(dx)^3)));

beta11=0.5*C11*(Q11^2+2*Q12^2)+C12*Q12*(2*Q11+Q12);
beta12=C11*Q12*(2*Q11+Q12)+C12*(Q11^2+3*Q12^2+2*Q11*Q12)+C44*Q44^2;
felas2=beta11*(P1(:,:,1:32).^4+P2(:,:,1:32).^4+P3(:,:,1:32).^4)+beta12*(P1(:,:,1:32).^2.*P2(:,:,1:32).^2+P3(:,:,1:32).^2.*P2(:,:,1:32).^2+P1(:,:,1:32).^2.*P3(:,:,1:32).^2);
Felas2=sum(sum(sum(felas2*(dx)^3)));

q11=C11*Q11+2*C12*Q12;
q12=C11*Q12+C12*(Q11+Q12);
q44=2*C44*Q44;
felas3=-(q11*yibu11+q12*yibu22+q12*yibu33).*P1(:,:,1:32).^2-(q11*yibu22+q12*yibu11+q12*yibu33).*P2(:,:,1:32).^2-...
    (q11*yibu33+q12*yibu11+q12*yibu22).*P3(:,:,1:32).^2-2*q44*(yibu12.*P1(:,:,1:32).*P2(:,:,1:32)+yibu23.*P2(:,:,1:32).*P3(:,:,1:32)+yibu13.*P2(:,:,1:32).*P3(:,:,1:32));
Felas3=sum(sum(sum(felas3*(dx)^3)));

felas=felas1+felas2+felas3;
felas=cat(3,felas,zeros(128,128,4));
Felas=sum(sum(sum(felas*(dx)^3)));

f=fbulk+fwall+felec+felas;

end

