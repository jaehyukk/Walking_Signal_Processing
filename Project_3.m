tic;
close all;clear all; clc;

path = './data';   
filename = 'eventdata.txt';

data = load(filename);
L = size(data,1);
SAMPLE_RATE = 32;

%% G단위 Scaling

g = 9.81; % G = 9.80665m/s^2
G_saturation = 41 * g; % G_saturation 범위 임의 설정
N_resolution = 10; % 양자화 범위 2^10-1 까지 표현
G_divisor = 2^(N_resolution + 1) / G_saturation; % G_divisor 계산
ac = (data./G_divisor); % data 값을 G_divisor로 나눈 데이터 ac에 저장
acc = (data./G_divisor)./100; 
% overflow 방지하기 위해 100으로 나눠서 데이터 acc에 저장

figure() 
hold on; grid on; box on;
plot(data); % data 그래프 플롯
xlabel('Number of Samples') % x축 지정
ylabel('G(9.8m/s^2)') % y축 지정
title('Raw Signal') % 제목 지정
legend('x','y','z') % 범례 지정
hold off; grid off; box off; 

figure() 
hold on; grid on; box on;
plot(ac); % ac 그래프 플롯
xlabel('Number of Samples') % x축 지정
ylabel('G(9.8m/s^2)') % y축 지정
title('scaling normalized signal') % 제목 지정
legend('x','y','z') % 범례 지정
hold off; grid off; box off; 

figure() 
hold on; grid on; box on;
plot(acc); % acc 그래프 플롯
xlabel('Number of Samples') % x축 지정
ylabel('G(9.8m/s^2)') % y축 지정
title('overflow prevention signal ') % 제목 지정
legend('x','y','z') % 범례 지정
hold off; grid off; box off;

%% 가속도 모멘텀 산출

for i = 1:L
momentum(i) = sqrt(acc(i,1).^2 + acc(i,2).^2 + acc(i,3).^2);
% momentum 산출, acc(i,1) = x축 데이터, acc(i,2) = y축 데이터, acc(i,3) = z축 데이터
end
m = mean(momentum); % 모멘텀의 평균값
result = momentum - m; % 결과 = 모멘텀 - 평균값

figure()
plot(result) % 모멘텀 결과값 그래프 플롯
legend('Momentum') % 범례 지정
xlabel('Number of Samples') % x축 지정
ylabel('G(9.8m/s^2)') % y축 지정
title('Momentum Signal in Samples') % 제목 지정

%% Low pass filter를 통한 노이즈

Sf = 32; % Sampling Frequency = 32Hz
n = 5; % Order = 5
cf = 0.5; % Cut-off Frequency = 0.5Hz

LPF = designfilt('lowpassfir', 'FilterOrder', n, 'CutoffFrequency', cf,'SampleRate', Sf);
% Low Pass Filter 적용

momentum_filter = filter(LPF,result); % LPF를 적용한 모멘텀 신호
figure()
plot(momentum_filter) % 노이즈 제거된 모멘텀 그래프 플롯
legend('FIR filter') % 범례 지정
xlabel('Number of Samples') % x축 지정
ylabel('G(9.8m/s^2)') % y축 지정
title('Momentum Signal in Samples') % 제목 지정

%% 모멘텀 신호 양자과정 (가우스함수)
Q_level = 0.005; % Q_level 지정
Q_bias = 1; % Q_bias 지정
m_quantization = floor((momentum_filter + Q_bias) / Q_level) * Q_level;
% 양자화 공식 적용

figure()
hold on; grid on; box on;
plot(m_quantization); % 양자화된 데이터 그래프 플롯
legend('Quantization') % 범례 지정
xlabel('Number of Samples') % x축 지정
ylabel('G(9.8m/s^2)') % y축 지정
title('Quantized Momentum Signal in Samples') % 제목 지정

%% Peak / Valley 검출

S = zeros(length(m_quantization),1); % S를 0으로 초기화
P = zeros(length(m_quantization),1); % P를 0으로 초기화

for i=1:length(m_quantization)-1 % S 구하기
    if m_quantization(i) < m_quantization(i+1) % 이전값 < 현재값
        S(i) = 1;
    elseif m_quantization(i) > m_quantization(i+1) % 이전값 > 현재값
        S(i) = 0;
    else % otherwise
        S(i+1) = S(i); % 이전값 대입
    end
end

for i = 2:length(S) % P 구하기
    if S(i-1) ~= S(i)
        if S(i-1) == 1 % S[i] > S[i-1]
            P(i) = 1; % Peak
        elseif S(i-1) == 0 % S[i] < S[i-1]
            P(i) = 0; % Valley
        end

        if i > 4 % 예외처리
        if (P(i) == 0) && (P(i-1) == 1)
            P(i) = NaN;
            P(i-1) = NaN;
        elseif (P(i) == 0) && (P(i-2) == 1)
            P(i) = NaN;
            P(i-2) = NaN;
        elseif (P(i) == 0) && (P(i-3) == 1)
            P(i) = NaN;
            P(i-3) = NaN;
        end
        end
    else
        P(i) = NaN; % 무효값 처리
    end
end

Peak = P; Valley = P; % peak, valley 길이 = 양자화된 데이터 길이

figure()
hold on; grid on; box on;
plot(find(P==1), m_quantization(find(P==1)),'ob') % Peak
plot(find(P==0), m_quantization(find(P==0)),'or') % Valley
ylim([0.9 1.04]) % y축 범위 조정
legend('Peak','Valley') % 범례 지정
xlabel('Number of Samples') % x축 이름 조정
ylabel('G(9.8m/s^2)') % y축 이름 조정
title('Peak/Valley Points in Samples') % 제목 조정

%% 걸음수 검출

Peak_value = m_quantization(find(P==1)); % peak 값 산출
Valley_value = m_quantization(find(P==0)); % valley 값 산출 

for i = 1:length(Peak_value) % peak와 valley 차이값 계산
    distance(i) = Peak_value(i) - Valley_value(i);
end

THRESHOLD_SWING = mean(distance)*0.5;  % THRESHOLD_SWING

P_index = find(P==1); % peak인 점의 index 대입
V_index = find(P==0); % valley인 점의 index 대입

THRESHOLD_ACTIVE = 1; % THRESHOLD_ACTIVE, 임의 조정

THRESHOLD_INTERVAL_TIME = zeros(length(P_index), 1);
for i = 1:length(P_index) % peak와 valley의 시간 간격
    THRESHOLD_INTERVAL_TIME(i) = P_index(i) - V_index(i);
end

THRESHOLD_INTERVAL_LONG_TIME = max(THRESHOLD_INTERVAL_TIME);
% peak와 valley 사이 시간 간격 최댓값
THRESHOLD_INTERVAL_SHORT_TIME = min(THRESHOLD_INTERVAL_TIME);
% peak와 valley 사이 시간 간격 최솟값
THRESHOLD_INTERVAL_LONG_SAMPLE = round(THRESHOLD_INTERVAL_LONG_TIME * SAMPLE_RATE);
% THRESHOLD_INTERVAL_LONG_SAMPLE
THRESHOLD_INTERVAL_SHORT_SAMPLE = round(THRESHOLD_INTERVAL_SHORT_TIME * SAMPLE_RATE);
% THRESHOLD_INTERVAL_SHORT_SAMPLE

i = 1; k =1;
for j = 1 : length(P_index) % j = 73까지 반복문 진행
    if i < length(P_index)
    i = i + 1; % i++
    if Peak_value(i) ~= -1 % 현재 peak ~= -1
        if Peak_value(i-1) >= THRESHOLD_ACTIVE % 이전 peak >= active
            if Peak_value(i) >= THRESHOLD_ACTIVE % 현재 peak >= active
                if Peak_value(i) - Valley_value(i-1) >= THRESHOLD_SWING
                    % 현재 peak - 이전 valley >= swing
                    if P_index(i) - P_index(i-1) <= THRESHOLD_INTERVAL_LONG_TIME
                        % 현재 peak - 이전 peak <= interval_long
                        if P_index(i) - P_index(i-1) >= THRESHOLD_INTERVAL_SHORT_TIME
                            % 현재 peak - 이전 peak <= interval_short
                            if round((P_index(i) - P_index(i-1)) * SAMPLE_RATE) <= THRESHOLD_INTERVAL_LONG_SAMPLE
                                % 현재 peak - 이전 peak <= interval_long
                                if round((P_index(i) - P_index(i-1)) * SAMPLE_RATE) >= THRESHOLD_INTERVAL_SHORT_SAMPLE
                                    % 현재 peak - 이전 peak <= interval_short
                                    Step(i) = 1;
                                else
                                    Step(i) = 0;
                                    no_Peak(k) = i;
                                    k = k + 1;
                                end
                            else
                                Step(i) = 0;
                                no_Peak(k) = i;
                                k = k + 1;
                            end
                        else
                            Step(i) = 0;
                            no_Peak(k) = i;
                            k = k + 1;
                        end
                    else
                        Step(i) = 0;
                        no_Peak(k) = i;
                        k = k + 1;
                    end
                else
                    Step(i) = 0;
                    no_Peak(k) = i;
                    k = k + 1;
                end
            else
                Step(i) = 0;
                no_Peak(k) = i;
                k = k + 1;
            end
        else
            Step(i) = 0;
            no_Peak(k) = i;
            k = k + 1;
        end
    else
        Step(i) = 0;
        no_Peak(k) = i;
        k = k + 1;
    end
  end
end

cnt=0;
for i = 1:length(P)
    if P(i) == 1
        cnt = cnt + 1;
    end
    if cnt == 2
        P(i) = 0;
    elseif cnt == 30
        P(i) = 0;
    elseif cnt == 31
        P(i) = 0;
    end
end

% 인식률 산출
N_r = 71;
N_a = sum(Step);
R = (N_r - (abs(N_r - N_a)))/N_r * 100;

disp(R); % 97.1831
disp(N_a); % 69
disp(THRESHOLD_ACTIVE); % 1
disp(THRESHOLD_SWING); % 0.0180
disp(THRESHOLD_INTERVAL_SHORT_TIME); % 12
disp(THRESHOLD_INTERVAL_LONG_TIME); % 34
disp(THRESHOLD_INTERVAL_SHORT_SAMPLE); % 384
disp(THRESHOLD_INTERVAL_LONG_SAMPLE); % 1088

figure()
hold on; grid on; box on;
plot(m_quantization)
plot(find(P==1), m_quantization(find(P==1)),'ob')
xlabel('Number of Samples')
ylabel('G(9.8m/s^2)')
title('WalkStep Detection in Samples')
legend('m quantization','Walk');

figure()
x = floor(P);
stem(x);
legend('Step');
toc;
