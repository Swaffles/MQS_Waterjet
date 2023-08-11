function myPSD = psdWaterJet(dataTable,testMatrix,testMatrixVars,dataFileName,time,debug)
% This function takes in a table containting raw data from a labview
% testing log called dataTable, the corresponding testMatrix, the 
% dataFileNames, as well as the option to set a debugger. The function
% queries the user to select the bounds for valid data. The function then 
% outputs a table called mypsd that contains the power and freq for Fx, Fy,
% Fz, Mx, My, Mz, Rotor Thrust, Rotor Speed (rps), and Net Force.
% E.g.
% Input dataTable:
% Fx    Fy      Fz      ...
% 101   2.99
% 99    3.03
% .     .       .       ...
% .     .       .       ...

% Output myTable:
% Fx Power      Fy Power    ...
% double        double      .
% double        double      .
% .             .           .
% .             .           .


maxfreq = 1000; % sampled at 2000 so 1/2
windowSz = (60/64)*2000; % 60 sec of data at 2000/sec, want ~ windows
nOverLap = 0.5*windowSz; % 10% overlap on window
% parse time
start0 = time(1);
last0 = time(2);
start = time(3);
last = time(4);

     %locate the trial in the test matrix incase missaligned
    testMatrixIndex = find(strcmp(dataFileName,testMatrix.TrialName));
    if debug
        fprintf("Building myFFT for %s\n", dataFileName);
        %fprintf("Using data from row %0.f of test matrix, corresponding to %s\n",...
        %        testMatrixIndex,testMatrix.TrialName(testMatrixIndex));
    end
    vars = fieldnames(dataTable);
    vars = vars(2:end-3); % drop last three table intrinsic properties
    vars(10) = []; % drop rotor-speed pulse count
    vars = strrep(vars,'_','-');
    vars = [vars(1:6);vars(8);vars(10)]; % the ones we want to keep for FFT-ing
    sz = length(vars);
    varNames = strings(sz,1);
    for i = 1:sz
        varNames(i) = strcat(vars(i)," Power");
    end
    varNames(end+1:end+2) = ["Total Net Force Power","Frequency"];
    varNames(end+1:end+2) = [testMatrixVars(1),testMatrixVars(2)];
    varTypes = ["double","double","double","double","double","double",...
        "double","double","double","double","double","double"];
    cz = length(varNames);
    % need to determine the min size for mypsd  
    tempPSD = pwelch(dataTable.(vars{1})(start:last),windowSz,nOverLap);
    %tempPSD(1) = [];
    n = length(tempPSD);
    tempPSD = 10*log10(tempPSD);
    height = length(tempPSD);
    myPSD = table('Size',[height,cz],'VariableTypes',varTypes,'VariableNames',varNames);
  
    %body forces and moments
    % get mean of when no jet is on
    ForceNoJet = mean(dataTable{start0:last0,[vars(1:6)]},1,"omitnan");
    
    BodyForce = dataTable{start:last,[vars(1:6)]};

    NetForce = BodyForce - ForceNoJet; % what will be written to Fx, Fy, ...
    
    TotalNetForce = sqrt(NetForce(:,1).^2+NetForce(:,2).^2+NetForce(:,3).^2);
    
    % compute FFTs
    [NetForcePSD,Freq] = pwelch(NetForce,hamming(windowSz),nOverLap,[],2*maxfreq);
    %NetForceFFT(1,:) =[]; % remove the first value (sum of the data)
    n = length(NetForcePSD);
    NetForce_Power = 10*log10(NetForcePSD);
    
    % frequency range
    %Freq = (1:n/2)/(n/2)*maxfreq;

    % total net force
    TotalNetForcePSD = pwelch(TotalNetForce,hamming(windowSz),nOverLap); %psd
    %TotalNetForcePSD(1) = [];
    n = length(TotalNetForcePSD);
    TotalNetForceFFT_Power = 10*log10(TotalNetForcePSD);

    % Rotor Reaction Force
    PSDRotorThrust = pwelch(dataTable.RotorThrust(start:last),hamming(windowSz),nOverLap);
    %PSDRotorThrust(1) = [];
    n = length(PSDRotorThrust);
    RotorThrust_Power = 10*log10(PSDRotorThrust);

    % Rotor speed
    PSDRotorSpeed = pwelch(dataTable.("RotorSpeed-Frequency")(start:last),hamming(windowSz),nOverLap);
    %FFTRotorSpeed(1) =[];
    n = length(PSDRotorSpeed);
    %RotorSpeed_Power = abs(FFTRotorSpeed(1:floor(n/2))).^2; 
    RotorSpeed_Power = 10*log10(PSDRotorSpeed);
    
    myPSD{:,varNames(1:6)} = NetForce_Power;
    myPSD{:,"RotorThrust Power"} = RotorThrust_Power;
    myPSD{:,"RotorSpeed-Frequency Power"} = RotorSpeed_Power;
    myPSD{:,"Total Net Force Power"} = TotalNetForceFFT_Power;
    myPSD{:,"Frequency"} = Freq;
    myPSD{:,"WaterjetSpeed"} = testMatrix.(testMatrixVars{1})(testMatrixIndex)*ones(floor(n),1);
    myPSD{:,"TrimAngle"} = testMatrix.(testMatrixVars{2})(testMatrixIndex)*ones(floor(n),1);

    %semilogy(Freq,RotorThrust_Power);
  
    clear temp
    return
end