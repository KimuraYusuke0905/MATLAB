%透磁率
m=4*pi*10^(-7);

%運転電流I[A]
I=110;

%メインコイルの諸元(入力内容)
z1m=10;%中心からの距離[mm]
N1m=104;%層数
N2m=28;%段数
a1m=70;%内半径[mm]

%メインコイルの諸元(自動出力)
Nm=N1m*N2m;%総巻き数
xm=0.2*N1m;%層長さ[mm]
ym=4*N2m;%段長さ[mm]
z2m=z1m+ym;%[mm]
a2m=a1m+xm;%外半径[mm]
Jm=(Nm*I)/(xm*ym*10^(-6));%電流密度[A/m]
am=(a2m)/(a1m);%α

%キャンセルコイルの諸元(入力内容)
z1c=175;%中心からの距離[mm]
N1c=90;%層数[mm]
N2c=6;%段数
a1c=70;%内半径[mm]

%キャンセルコイルの諸元(自動出力)
Nc=N1c*N2c;%総巻き数
xc=0.2*N1c;%層長さ[mm]
yc=4*N2c;%段長さ[mm]
z2c=z1c+yc;%[mm]
a2c=a1c+xc;%外半径[mm]
Jc=(Nc*I)/(xc*yc*10^(-6));%電流密度[A/m]
ac=(a2c)/(a1c);%α


%位置z [mm]
t=linspace(-200,200,401);%0.4[mm]刻み
z=transpose(t);

%メインコイルの磁場計算
b1m=(z1m-z)/a1m;%β1
b2m=(z2m-z)/a1m;%β2
F1m=m.*b1m.*(log((am+sqrt(am.^2+b1m.^2))./(1+sqrt(1+b1m.^2))));
F2m=m.*b2m.*(log((am+sqrt(am.^2+b2m.^2))./(1+sqrt(1+b2m.^2))));
Fm=F2m-F1m;
X1=Fm(z<=0);
X2=flip(Fm(z>=0));
Xa=X1+X2;
Xb=flip(Xa);
Xb(1,:)=[];
X=[Xa;Xb];
Bm=Jm.*a1m*10^(-3)*X;

%キャンセルコイルの磁場計算
b1c=(z1c-z)/a1c;%β1
b2c=(z2c-z)/a1c;%β2
F1c=m.*b1c.*(log((ac+sqrt(ac.^2+b1c.^2))./(1+sqrt(1+b1c.^2))));
F2c=m.*b2c.*(log((ac+sqrt(ac.^2+b2c.^2))./(1+sqrt(1+b2c.^2))));
Fc=F2c-F1c;
Y1=Fc(z<=0);
Y2=flip(Fc(z>=0));
Ya=Y1+Y2;
Yb=flip(Ya);
Yb(1,:)=[];
Y=[Ya;Yb];
Bc=Jc.*a1c*10^(-3)*Y;

%全体の磁場計算
B=Bm-Bc;

%dB/dz
dBdz=gradient(B,z);

%グラフへの描写
set(0, 'DefaultLineLineWidth', 1);
yyaxis left;
plot(z,Bm,'r--',z,Bc,'b--',z,B,'k-');
xlabel('位置z [mm]');
ylabel('磁場B [T]');
title('Magnetic field calculation');
yyaxis right;
plot(z,dBdz,'c-');
ylabel('dB/dz [T/mm]');
legend('Main','Shield','Mix');





