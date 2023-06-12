clear;
clc;
close all;
%% Basic Electromagnetic Parameters
Frequency = 10e9;
Lightspeed = physconst('LightSpeed');%光速
Wavelength = Lightspeed/Frequency;%波长
Wavenumber = 2*pi/Wavelength;%
%% Array Parameters
N =33;
X = (-(N-1)/2:1:(N-1)/2)*Wavelength/2;
%X(5) = 0*X(5);
I =  ones(1,N);
%n = N;
%A = [0,1];
%I = A(randi(numel(A),1,n));
theta0 = 0;
alpha = -Wavenumber*X*sind(theta0);
G = exp(1j*alpha);
%我需要把alpha做成对于0对称的情况:
%% ArrayFactor Samping
Ns =500;% Sampling number
theta = linspace(-90,90,Ns);
E =zeros(1,Ns);
win=chebwin(N);%加窗
%%fomat
%随机化处理，将I1进行调整
bw = chebwin(N);
I1 = I;
I1(10) = -1;
I1(2) = -1;
%% 新的调质处理
I2 = I1;
I2(10)=0;
I2(2)=0;
%I2 = I2.*bw';
%% 新的处理
I3 = I1;
I3(20)=-1;
I3(22)=-1;
%% 新的处理
I4 = I2;
I4(20)=0;
I4(22)=0;
%% add noise：
SNR = 10;
INR = 10;
noise=[randn(N,Ns)+1i*randn(N,Ns)]/sqrt(2);
noise=noise';%噪声成功添加，计算:
thetaj=[20,40,70];%噪声信号方向(零陷方向)
nj = length(thetaj);
as0=exp(Wavenumber*j*pi*[0:N-1]'*sind(theta0));
amp_j=db(INR)*0.707*(randn(nj,Ns)+1i*randn(nj,Ns));
%amp_j1 = db(INR)*0.707*(randn(nj,Ns)+1i*randn(nj,Ns));
aj=exp(j*pi*[0:N-1]'*sin(thetaj*pi/180));
%%这里区分出只使用as0的信号用以对比
%xj1=as0*amp_j+noise';
xj=aj*amp_j+noise';
Rin=1/Ns*(xj*xj');
for m= 1:1:N
    for n=1:1:N
        Rin(m,n) = Rin(m,n)*sinc((m-n)*5)*pi/3;
    end
end
Rin_inv=inv(Rin);
W=Rin_inv*as0;
%%得再做一个零陷扩展:
for num = 1:Ns
    E(num)=sum(I.*exp(1j*(Wavenumber*X*(sind(theta(num))-sind(theta0))))+1e-3;
    E(num)=sum(I.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*
