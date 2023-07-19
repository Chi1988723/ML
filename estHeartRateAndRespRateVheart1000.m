% Program for HeartRate & RespirationRate estimation
%--------------------------------------------------------------------------
% input

close all;
clear;

% column 1: time; column 2: 24 GHz radar I-channel; column 3: 24 GHz radar Q-channel; column 4: 10 GHz radar I-channel; column 5: respiratory band signal; and column 6: ECG signal.

fs      = 1000 ;     % Sampling rate [Hz]
start   = 0 ;     % Seconds to start data [s]
stop    = 600 ;     % Seconds to end data [s]
desample = 50;
deFs = fs / desample;
lenSection = 10;
nSection = 60;

lpFilt = designfilt('highpassfir', 'PassbandFrequency', 1.8 / deFs * 2, ...% 0.4 * 
                  'StopbandFrequency', 1.2 / deFs * 2, 'PassbandRipple', 3, ...
                  'StopbandAttenuation', 30);

MSE = zeros(1,9);
result = zeros(2, nSection, 9);
for iii = 1 : 9
    file = readmatrix(['subject',num2str(iii),'.csv']) ;
    Data = file(22:end,1:6) ;
    signalI = Data(start*fs+1:stop*fs, 2);
    signalI = signalI(1 : desample : end);
    signalQ = Data(start*fs+1:stop*fs,3) ;
    signalQ = signalQ(1 : desample : end);
    radarSig = Data(start*fs+1:stop*fs,6);
    ecgSig = radarSig(1 : 1 : end);
    signalC = signalI + 1j * signalQ;

    signalI = filter(lpFilt, signalC);

    signalHandle.isPlot = 0;
    signalHandle.fftLenMuti = 10;
    signalHandle.nInt = 10;
    defaultFreRange = [0.75; 1.5];
    signalHandle.periodRange = flipud( 1 ./ defaultFreRange);
    for ii = 1 : nSection       
       % 参考的ECG信号
       
       signalHandle.fs = fs;
       signalHandle.signal = ecgSig(lenSection * ((ii - 1) * fs) + 1 : lenSection * (ii * fs));
       signalHandle.freRange = defaultFreRange;
       referFre = myFFT(signalHandle);
       result(1, ii, iii) = 60 * referFre;

       signalHandle.fs = deFs;
       signalHandle.signal =  signalI(lenSection * ((ii - 1) * deFs) + 1 : lenSection * (ii * deFs));
       temp = myPeriodML(signalHandle);
       MSE(1, iii) = MSE(1, iii) + (referFre - 1 / temp) .^ 2; 
       result(2, ii, iii) = 60 / temp;


    end

%     figure(1); clf;
%     plot(result( : , : , iii)');
end
legend('PAMDF', 'AMDF', 'Best match');

MSE = 60 * sqrt(MSE / nSection);
figure;
boxplot(MSE');
mean(MSE, 2)





