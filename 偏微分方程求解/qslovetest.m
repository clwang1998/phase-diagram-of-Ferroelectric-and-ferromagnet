C11=1.78*10^11;
C12=0.964*10^11;
C44=1.22*10^11;
a11=ones(128,128);
a12=a11;
a13=a11;
a21=a11;
a22=a11;
a23=a11;
a31=a11;
a32=a11;
a33=a11;
a41=a11;
a42=a11;
a43=a11;
a51=a11;
a52=a11;
a53=a11;
a61=a11;
a62=a11;
a63=a11;
D1=a11;
D2=a11;
D3=a11;
D4=a11;
D5=a11;
D6=a11;
bian1=ones(128,128);
bian2=ones(128,128);
bian3=ones(128,128,32);
hf=1;
yinta=sqrt(bian1.^2+bian2.^2);
cu11f=a11.*exp(-yinta*hf);
cu12f=a12.*exp(-yinta*hf);
cu13f=a13.*exp(-yinta*hf);
cu21f=a21.*exp(-yinta*hf);
cu22f=a22.*exp(-yinta*hf);
cu23f=a23.*exp(-yinta*hf);
cu31f=a31.*exp(-yinta*hf)+1i*hf*a21.*exp(-yinta*hf);
cu32f=a32.*exp(-yinta*hf)+1i*hf*a22.*exp(-yinta*hf);
cu33f=a33.*exp(-yinta*hf)+1i*hf*a23.*exp(-yinta*hf);
cu41f=a41.*exp(yinta*hf);
cu42f=a42.*exp(yinta*hf);
cu43f=a43.*exp(yinta*hf);
cu51f=a51.*exp(yinta*hf);
cu52f=a52.*exp(yinta*hf);
cu53f=a53.*exp(yinta*hf);
cu61f=a61.*exp(yinta*hf)+1i*hf*a51.*exp(yinta*hf);
cu62f=a62.*exp(yinta*hf)+1i*hf*a52.*exp(yinta*hf);
cu63f=a63.*exp(yinta*hf)+1i*hf*a53.*exp(yinta*hf);
eu11=-yinta.*a11.*exp(-yinta*hf);
eu12=-yinta.*a12.*exp(-yinta*hf);
eu13=-yinta.*a13.*exp(-yinta*hf);
eu21=-yinta.*a21.*exp(-yinta*hf);
eu22=-yinta.*a22.*exp(-yinta*hf);
eu23=-yinta.*a23.*exp(-yinta*hf);
eu31=-yinta.*a31.*exp(-yinta*hf)+1i*a21.*exp(-yinta*hf)-1i*hf*yinta.*a21.*exp(-yinta*hf);
eu32=-yinta.*a32.*exp(-yinta*hf)+1i*a22.*exp(-yinta*hf)-1i*hf*yinta.*a22.*exp(-yinta*hf);
eu33=-yinta.*a33.*exp(-yinta*hf)+1i*a23.*exp(-yinta*hf)-1i*hf*yinta.*a23.*exp(-yinta*hf);
eu41=yinta.*a41.*exp(yinta*hf);
eu42=yinta.*a42.*exp(yinta*hf);
eu43=yinta.*a43.*exp(yinta*hf);
eu51=yinta.*a51.*exp(yinta*hf);
eu52=yinta.*a52.*exp(yinta*hf);
eu53=yinta.*a53.*exp(yinta*hf);
eu61=yinta.*a51.*exp(yinta*hf)+1i*a51.*exp(yinta*hf)+1i*hf*yinta.*a51.*exp(yinta*hf);
eu62=yinta.*a52.*exp(yinta*hf)+1i*a52.*exp(yinta*hf)+1i*hf*yinta.*a52.*exp(yinta*hf);
eu63=yinta.*a53.*exp(yinta*hf)+1i*a53.*exp(yinta*hf)+1i*hf*yinta.*a53.*exp(yinta*hf);
hs=1;
hss=-hs;
cu11s=a11.*exp(-yinta*hss);
cu12s=a12.*exp(-yinta*hss);
cu13s=a13.*exp(-yinta*hss);
cu21s=a21.*exp(-yinta*hss);
cu22s=a22.*exp(-yinta*hss);
cu23s=a23.*exp(-yinta*hss);
cu31s=a31.*exp(-yinta*hss)+1i*hss*a21.*exp(-yinta*hss);
cu32s=a32.*exp(-yinta*hss)+1i*hss*a22.*exp(-yinta*hss);
cu33s=a33.*exp(-yinta*hss)+1i*hss*a23.*exp(-yinta*hss);
cu41s=a41.*exp(yinta*hss);
cu42s=a42.*exp(yinta*hss);
cu43s=a43.*exp(yinta*hss);
cu51s=a51.*exp(yinta*hss);
cu52s=a52.*exp(yinta*hss);
cu53s=a53.*exp(yinta*hss);
cu61s=a61.*exp(yinta*hss)+1i*hss*a51.*exp(yinta*hss);
cu62s=a62.*exp(yinta*hss)+1i*hss*a52.*exp(yinta*hss);
cu63s=a63.*exp(yinta*hss)+1i*hss*a53.*exp(yinta*hss);
aa=1;
while aa<129
    bb=1;
    while bb<129
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
bb=bb+1;
    end
    aa=aa+1;
end
q1=Q(:,:,1);
q2=Q(:,:,2);
q3=Q(:,:,3);
q4=Q(:,:,4);
q5=Q(:,:,5);
q6=Q(:,:,6);

