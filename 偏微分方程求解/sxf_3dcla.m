function sxf_3dcla
clear;clc;
for z=40
    z1=10*z;
    str1='F:/Datas/Single/finite/3D/PW/4nm';
    str2='/particles';
    str3=num2str(z,'%04d');
    str4='.sdf';
    str=strcat(str1,str2,str3,str4);
    [b,h]=lv(str);
    x=-18.6/0.8:50/1000/0.8:18.6/0.8-8/1000/0.8;
    y=-2.4/0.8:16/1000/0.8:10.4/0.8;
    z=-18.6/0.8:50/1000/0.8:18.6/0.8-8/1000/0.8;
for i=91
choose=i;
left=1;
right=400;
switch choose
    case 1
        ex=gd(b,h,'ex')/4.02*10^-12;
        ex1=ex(:,:,372);
        figure
        imagesc(ex1')
        colorbar
        ex2=ex(200:1800,:,:);
        [x,y,z]=meshgrid(x,y,z);
        p=patch(isosurface(x,y,z,ex2,5));view(105,30);
        set(p,'FaceColor','red','EdgeColor','none');
        camlight
        lighting gouraud
        
        
        np21=np(:,:,372)*100+8;
        np21(np21<13)=0;
        ne21=ne20*10+32;
        ne21(ne21<35)=0;
        all=ex21+np21+ne21;
        zhi=-12*ones(size(xx));
        hold on;
        surf(xx,yy,zhi,all,'edgecolor','none','facecolor','interp')
        
    case 91
        tic
        %figure
        r=6*10^-6;
        aa='proton';ab='carb';ac='electron1';ad='gold';strr='subset_bw/';
        strr0=strcat(strr,aa);
        strr1=strcat('gamma/',strr0);strr2=strcat('grid/',strr0);strr3=strcat('weight/',strr0);
        strr4=strcat('px/',strr0);strr5=strcat('py/',strr0);
        gam=gd(b,h,strr1);
        ygam=gd(b,h,strr2);ygam1=ygam.y;xgam1=ygam.x;zgam1=ygam.z;rgam1=sqrt(ygam1.^2+zgam1.^2);
        ppx=gd(b,h,strr4);ppy=gd(b,h,strr5);
        theta=ppy./ppx;theta1=theta*180/pi;
        ind=find(rgam1>-r&rgam1<r&abs(ygam1)>-.1*10^-6&abs(theta1)<90);%&xgam1>-100e-6&xgam1<250e-6);
        gamm=gam(ind);
        wei=gd(b,h,strr3);
        weim=wei(ind);
        Ep=(gamm-1)*939.375;
        dE=.5;count=0.0;
        En=[];dN1=[];
        for k=1:1000
            in=find(Ep>=k*dE&Ep<(k+1)*dE);
            weim1=weim(in);
            %weight2=weight1;
            dN=sum(weim1);
            dN1(k)=dN;
            En(k)=k*dE;
        end
        %plot(En,dN1,'.');
        dN2=dN1*pi*r/2*10/dE;
        semilogy(En,dN1*20/dE,'m','LineWidth',2);hold on
        xlabel('En [MeV/u]','fontsize',16);ylabel('dN/dE','fontsize',16);
        %xlim([0 100]);
        %ylim([10^12 10^13]);
        title(strcat('Energy Spectra of Proton (t=',num2str(z1),'T0)'),'fontsize',16);
        %print('-dpng',['F:/Datas/Test/RPA/square/heavy/5um/pl/4/Figures/En\Enc_' num2str(z1) '.png']);  
        toc
        
        
end
end
end
end