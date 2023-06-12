clear;
clc;
close all;
%% Basic Electromagnetic Parameters
Frequency = 10e9;
Lightspeed = physconst('LightSpeed');%光速
Wavelength = Lightspeed/Frequency;%波长
Wavenumber = 2*pi/Wavelength;%
%% Array Parameters
N =101;
X = (-(N-1)/2:1:(N-1)/2)*Wavelength/2;
%X(5) = 0*X(5);
I =  ones(1,N);
%n = N;
%A = [0,1];
%I = A(randi(numel(A),1,n));
theta0 = [0];
theta1 = [15];
alpha = -Wavenumber*X*sind((theta0));
alpha1 = -Wavenumber*X*sind((theta1));
G = exp(1j*alpha);
%我需要把alpha做成对于0对称的情况:
%% ArrayFactor Samping
Ns =500;% Sampling number
theta = linspace(-90,90,Ns);
E =zeros(1,Ns);
win=blackman(N);%加窗,加窗效果现在不明显，头疼
%%fomat
%随机化处理，将I1进行调整
bw = chebwin(N);

I1 = I;
I1(10) = -1;
I1(2) = -1;
I1(40) = -1;
I1(42) = -1;

%% 新的调质处理
I2 = I1;
I2(10) = -0;
I2(2) = -0;
I2(40) = -0;
I2(42) = -0;
I1 = I1.*bw;
I2 = I2.*bw;
I = I;
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
thetaj=[20,40,60];%噪声信号方向(零陷方向)
nj = length(thetaj);
as0=exp(Wavenumber*j*pi*[0:N-1]'*sind(theta0(1)));
as1=exp(Wavenumber*j*pi*[0:N-1]'*sind(theta1(1)));
amp_j=db(INR)*0.707*(randn(nj,Ns)+1i*randn(nj,Ns));
%amp_j1 = db(INR)*0.707*(randn(nj,Ns)+1i*randn(nj,Ns));
aj=exp(j*pi*[0:N-1]'*sin(thetaj*pi/180));
%%这里区分出只使用as0的信号用以对比
%xj1=as0*amp_j+noise';
xj=aj*amp_j+noise';
Rin=1/Ns*(xj*xj');
for m= 1:1:N
    for n=1:1:N
        Rin(m,n) = Rin(m,n)*sinc((m-n)*15)*pi/3;
    end
end
Rin_inv=inv(Rin);
W=Rin_inv*as0;
%%得再做一个零陷扩展:
for num = 1:Ns
    %E(num)=sum(I.*exp(1j*(Wavenumber*X*(sind(theta(num))-sind(theta0))*win)))+1e-3;
    E(num)=sum(W'.*I.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*as0)+1e-3;
    E8(num)=sum(I.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*as1)+1e-3;
end
for num = 1:Ns
    %这里原本代码有一个加窗处理，存疑保留于此:
    %E1(num)=sum(W'.*I1.*(exp(1j*(Wavenumber*X*(sind(theta(num)))))))+1e-3;
    E1(num)=sum(W'.*I1.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*as0)+1e-3;
    E2(num)=sum(W'.*I2.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*as0)+1e-3;
    E3(num)=sum(W'.*I3.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*as0)+1e-3;
    E4(num)=sum(W'.*I4.*(exp(1j*(Wavenumber*X*(sind(theta(num))))))*as0)+1e-3;
end
E_dB = db(E)-max(db(E));
E_dB8 = db(E8)-max(db(E8));
E_db1 = db(E1)-max(db(E1));
E_db2 = db(E2)-max(db(E2));
E_db3 = db(E3)-max(db(E3));
E_db4 = db(E4)-max(db(E4));
%% plot figure
figure()
plot(theta,E_db1);
hold on;
plot(theta,E_db2);
hold on;
%plot(theta,E_db3);
%hold on;
%plot(theta,E_db4);
hold on;
plot(theta,E_dB);%normalized
hold on;
%plot(theta,E_dB8);%normalized
axis tight
grid on 
xlabel('角度/degree');ylabel('波束图/db');
legend("反转天线子集","开关天线子集","主瓣方向波束成形");
figure()
% Prepare data for QPSK modulation
data = randi([0, 1], 1, Ns); % Random data
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

% QPSK modulation
modulatedData = qpskMod(data.');

% Channel response for E_dB, E_db1, and E_db2
channelResponse_dB = E_dB(1:2:end) + 1j*E_dB(2:2:end);
channelResponse_db1 = E_db1(1:2:end) + 1j*E_db1(2:2:end);
channelResponse_db2 = E_db2(1:2:end) + 1j*E_db2(2:2:end);

% Apply channel response
received_dB = modulatedData .* channelResponse_dB.';
received_db1 = modulatedData .* channelResponse_db1.';
received_db2 = modulatedData .* channelResponse_db2.';

% Plot Constellations
figure()
plot(received_dB, 'o', 'DisplayName', '传统波束成形技术');
hold on;
plot(received_db1, 'o', 'DisplayName', '反转天线波束成形');
hold on;
plot(received_db2, 'o', 'DisplayName', '开关天线波束成形');
title('QPSK Constellations');
xlabel('In-phase'); ylabel('Quadrature');
legend('show');
grid on;
axis equal;

% 误码率计算参数
%% Prepare data for QPSK modulation
data = randi([0, 1], 1, Ns); % Random data
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

%% QPSK modulation
modulatedData = qpskMod(data.');

%% Prepare arrays for BER vs. angle plot
BER_dB_vs_angle = zeros(1, Ns);
BER_db1_vs_angle = zeros(1, Ns);
BER_db2_vs_angle = zeros(1, Ns);

%% Loop through angles and calculate BER
for i = 1:Ns
    % Channel response for E_dB, E_db1, and E_db2
    channelResponse_dB = E_dB(i) + 1j * E_dB(i);
    channelResponse_db1 = E_db1(i) + 1j * E_db1(i);
    channelResponse_db2 = E_db2(i) + 1j * E_db2(i);

    % Apply channel response
    received_dB = modulatedData .* channelResponse_dB;
    received_db1 = modulatedData .* channelResponse_db1;
    received_db2 = modulatedData .* channelResponse_db2;

    % QPSK demodulation and bit error rate calculation
    demodulatedData_dB = qpskDemod(received_dB);
    demodulatedData_db1 = qpskDemod(received_db1);
    demodulatedData_db2 = qpskDemod(received_db2);

    % Calculate bit error rate (BER)
    numErrors_dB = sum(demodulatedData_dB ~= data.');
    numErrors_db1 = sum(demodulatedData_db1 ~= data.');
    numErrors_db2 = sum(demodulatedData_db2 ~= data.');

    BER_dB_vs_angle(i) = numErrors_dB / Ns;
    BER_db1_vs_angle(i) = numErrors_db1 / Ns;
    BER_db2_vs_angle(i) = numErrors_db2 / Ns;
end

%% Plot BER vs. angle
figure();
plot(theta, BER_dB_vs_angle, 'LineWidth', 2, 'DisplayName', '传统波束成形');
hold on;
plot(theta, BER_db1_vs_angle, 'LineWidth', 2, 'DisplayName', '反转天线波束成形');
hold on;
plot(theta, BER_db2_vs_angle, 'LineWidth', 2, 'DisplayName', '开关天线波束成形');
xlabel('Angle (degree)');
ylabel('Bit Error Rate');
title('BER vs. Angle for QPSK modulation');
grid on;
legend('show');
%% Prepare data for QPSK modulation
data = randi([0, 1], 1, Ns); % Random data
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

%% QPSK modulation
modulatedData = qpskMod(data.');

%% Calculate channel capacity for E_dB
noiseVar_dB = 10.^(-SNR/10);
capacity_dB = log2(1 + abs(E_dB).^2 .* noiseVar_dB / 2);
rate_dB = capacity_dB / log2(4); % Log2(4) is the number of bits per symbol for QPSK modulation

%% Calculate channel capacity for E_db1
noiseVar_db1 = 10.^(-SNR/10);
capacity_db1 = log2(1 + abs(E_db1).^2 .* noiseVar_db1 / 2);
rate_db1 = capacity_db1 / log2(4); % Log2(4) is the number of bits per symbol for QPSK modulation

%% Calculate channel capacity for E_db2
noiseVar_db2 = 10.^(-SNR/10);
capacity_db2 = log2(1 + abs(E_db2).^2 .* noiseVar_db2 / 2);
rate_db2 = capacity_db2 / log2(4); % Log2(4) is the number of bits per symbol for QPSK modulation

%% Plot channel capacity vs. angle
figure();
plot(theta, rate_dB, 'LineWidth', 2, 'DisplayName', '反转波束形成安全速率');
hold on;
plot(theta, rate_db1, 'LineWidth', 2, 'DisplayName', '开关天线子集安全速率');
hold on;
plot(theta, rate_db2, 'LineWidth', 2, 'DisplayName', '传统波束成形安全速率');
xlabel('Angle (degree)');
ylabel('Capacity (bps/Hz)');
title('Secrecy rate vs. Angle for QPSK modulation');
grid on;
legend('show');
%% Prepare data for QPSK modulation
% Initialize variables
numBits = 1e5; % number of bits
theta = -90:2:90; % angles
Ns = length(theta); % number of angles
EbNo_dB = 10; % Eb/No in dB

% Generate random binary data
data = randi([0 1], numBits, 1);

% Create QPSK modulator and demodulator objects
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput', true);

% Prepare arrays for BER vs. angle plot
BER_dB_vs_angle = zeros(1, Ns);
BER_db1_vs_angle = zeros(1, Ns);
BER_db2_vs_angle = zeros(1, Ns);

% Loop through angles and calculate BER
for i = 1:Ns
    % Channel response for E_dB, E_db1, and E_db2
    channelResponse_dB = E_dB(i);
    channelResponse_db1 = E_db1(i);
    channelResponse_db2 = E_db2(i);

    % Apply channel response
    transmittedData = qpskMod(data);
    received_dB = awgn(transmittedData .* channelResponse_dB, EbNo_dB, 'measured');
    received_db1 = awgn(transmittedData .* channelResponse_db1, EbNo_dB, 'measured');
    received_db2 = awgn(transmittedData .* channelResponse_db2, EbNo_dB, 'measured');

    % QPSK demodulation and bit error rate calculation
    demodulatedData_dB = qpskDemod(received_dB);
    demodulatedData_db1 = qpskDemod(received_db1);
    demodulatedData_db2 = qpskDemod(received_db2);

    % Calculate bit error rate (BER)
    numErrors_dB = sum(demodulatedData_dB ~= data);
    numErrors_db1 = sum(demodulatedData_db1 ~= data);
    numErrors_db2 = sum(demodulatedData_db2 ~= data);

    BER_dB_vs_angle(i) = numErrors_dB / numBits;
    BER_db1_vs_angle(i) = numErrors_db1 / numBits;
    BER_db2_vs_angle(i) = numErrors_db2 / numBits;
end

% Plot BER vs. angle
figure();
plot(theta, BER_dB_vs_angle, 'LineWidth', 2, 'DisplayName', 'E_dB');
hold on;
plot(theta, BER_db1_vs_angle, 'LineWidth', 2, 'DisplayName', 'E_db1');
hold on;
plot(theta, BER_db2_vs_angle, 'LineWidth', 2, 'DisplayName', 'E_db2');
xlabel('Angle (degree)');
ylabel('Bit Error Rate');
title(sprintf('BER vs. Angle for QPSK modulation (Eb/No = %d dB)', EbNo_dB));
grid on;
legend('show');

figure();
%% Simulate BER for QPSK modulation

% 从之前的代码中获取已经计算好的波束方向图
% E_db1、E_db2 和 E_dB

% 仿真参数
EbN0_dB = -5:1:20; % 信噪比范围 (dB)
numSymbols = 1e6; % 符号数
QPSK = comm.QPSKModulator('BitInput', true); % QPSK调制器
QPSK_Demod = comm.QPSKDemodulator('BitOutput', true); % QPSK解调器
BER = zeros(3, length(EbN0_dB)); % 初始化误码率矩阵

% 生成随机比特
bits = randi([0 1], 1, 2*numSymbols);

% 蒙特卡罗仿真
for j = 1:length(EbN0_dB)
    fprintf('EbN0 = %d dB\n', EbN0_dB(j));
    % QPSK调制
    symbols = step(QPSK, bits');

    % 添加加权噪声
    noisePower = 10.^(-EbN0_dB(j)/10);
    noise = sqrt(noisePower/2) * (randn(size(symbols)) + 1i*randn(size(symbols)));

    % 为三个波束成形方向图分别计算误码率
    % 修改的代码段
% 修改的代码段
for k = 1:3
    if k == 1
        beamPattern = E_db1;
    elseif k == 2
        beamPattern = E_db2;
    else
        beamPattern = E_dB;
    end

    % 将波束方向图转换为复数权重向量
    beamPattern = 10.^(beamPattern/20);

    % 对波束权重进行插值，以使其与信号长度相匹配
    beamPattern_interp = interp1(1:length(beamPattern), beamPattern, linspace(1, length(beamPattern), length(symbols)), 'linear');

    % 对信号进行波束成形处理
    beamformedSymbols = symbols .* beamPattern_interp(:);

    % 添加噪声
    noisySymbols = beamformedSymbols + noise;

    % QPSK解调
    demodBits = step(QPSK_Demod, noisySymbols);

    % 计算误码率
    [~, ber] = biterr(bits', demodBits);
    BER(k, j) = ber;
end

end

% 画图
figure;
semilogy(EbN0_dB, BER(1, :), '-o', 'LineWidth', 2, 'DisplayName', '反转天线子集');
hold on;
semilogy(EbN0_dB, BER(2, :), '-x', 'LineWidth', 2, 'DisplayName', '开关天线子集');
semilogy(EbN0_dB, BER(3, :), '-s', 'LineWidth', 2, 'DisplayName', '主瓣方向波束成形');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('误码率（BER）');
legend('Location', 'southwest');
title('QPSK调制的误码率');
% 从之前的代码中获取已经计算好的波束方向图
% E_db1、E_db2 和 E_dB

% 仿真参数
% 从之前的代码中获取已经计算好的波束方向图
% E_db1、E_db2 和 E_dB

% 仿真参数
theta_deg = -90:3:90; % 角度范围 (度)
EbN0_dB = 13; % 信噪比 (dB) 调整为 0 dB
numSymbols = 1e5; % 符号数

QPSK = comm.QPSKModulator('BitInput', true); % QPSK调制器
QPSK_Demod = comm.QPSKDemodulator('BitOutput', true); % QPSK解调器
BER = zeros(3, length(theta_deg)); % 初始化误码率矩阵

% 生成随机比特
bits = randi([0 1], 1, 2*numSymbols);

% 为每个角度进行蒙特卡罗仿真
for i = 1:length(theta_deg)
    fprintf('角度 = %d 度\n', theta_deg(i));

    % QPSK调制
    symbols = step(QPSK, bits');

    % 添加加权噪声
    noisePower = 10.^(-EbN0_dB/10);
    noise = sqrt(noisePower/2) * (randn(size(symbols)) + 1i*randn(size(symbols)));

    % 为三个波束成形方向图分别计算误码率
    for k = 1:3
        if k == 1
            beamPattern = E_db1;
        elseif k == 2
            beamPattern = E_db2;
        else
            beamPattern = E_dB;
        end

        % 将波束方向图转换为复数权重向量
        beamPattern = 10.^(beamPattern/20);

        % 获取当前角度的波束权重
        beam_weight = interp1(linspace(-90, 90, length(beamPattern)), beamPattern, theta_deg(i), 'linear');

        % 对信号进行波束成形处理
        beamformedSymbols = symbols .* beam_weight;

        % 添加噪声
        noisySymbols = beamformedSymbols + noise;

        % QPSK解调
        demodBits = step(QPSK_Demod, noisySymbols);

        % 计算误码率
        [~, ber] = biterr(bits', demodBits);
        BER(k, i) = ber;
    end
end

% 画图
figure;
semilogy(theta_deg, BER(1, :), '-o', 'LineWidth', 2, 'DisplayName', '反转天线子集');
hold on;
semilogy(theta_deg, BER(2, :), '-x', 'LineWidth', 2, 'DisplayName', '开关天线子集');
semilogy(theta_deg, BER(3, :), '-s', 'LineWidth', 2, 'DisplayName', '主瓣方向波束成形');
grid on;
xlabel('角度 (度)');
ylabel('误码率（BER）');
legend('Location', 'southwest');
title('误码率随角度变化');






