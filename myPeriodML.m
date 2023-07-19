function [pamdf] = myPAMDF(signalHandle)

% 


    signal = signalHandle.signal;
    nInt = signalHandle.nInt;
    fs = signalHandle.fs;
    periodRange = signalHandle.periodRange;
    

    signal = signal - mean(signal); 
    signalI = real(signal);
    signalQ = imag(signal); 
    
    %% 

    delayScaleAMDF = (0 : length(signal) - 1)';
    amdf = zeros(length(delayScaleAMDF), 1);
    for ii = 1 : length(delayScaleAMDF)
        temp1 = abs(signalI(1 : end - delayScaleAMDF(ii)) - signalI(1 + delayScaleAMDF(ii) : end)) .^ 2;
        temp2 = abs(signalQ(1 : end - delayScaleAMDF(ii)) - signalQ(1 + delayScaleAMDF(ii) : end)) .^ 2;
        amdf(ii) = sum(temp1) + sum(temp2);
    end
    delayScaleFull = ((periodRange(1) * fs) : 1 / nInt : (periodRange(2) * fs)).';


    %% 
    nAMDF = length(amdf);
    fullAMDF = zeros(length(delayScaleFull), 1);
    for ii = 1 : length(delayScaleFull)
        weight = floor((nAMDF - 1) / delayScaleFull(ii)) * (nAMDF + 1) - sum(delayScaleAMDF(round(1 + delayScaleFull(ii) : delayScaleFull(ii) : nAMDF))) + 1;
%         weight = (nAMDF / delayScaleFull(ii) - 1);
        fullAMDF(ii) = sum(amdf(round(1 + delayScaleFull(ii) : delayScaleFull(ii) : nAMDF))) / weight;
    end
    
    
    [~, resultFre] = min(fullAMDF);
    pamdf = delayScaleFull(resultFre) / fs;
    
    %%
%     figure(2);
%     plot(delayScaleFull, fullAMDF);
 
end

