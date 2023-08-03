function myFFT = fftWaterJet(dataTable,testMatrix,testMatrixVars,dataFileName,time,debug)
% This function takes in a table containting raw data from a labview
% testing log called dataTable, the corresponding testMatrix, the 
% dataFileNames, as well as the option to set a debugger. The function
% queries the user to select the bounds for valid data. The function then 
% outputs a table called myFFT that contains the power and freq for Fx, Fy,
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
    % need to determine the min size for myFFT  
    tempFFT = fft(dataTable.(vars{1})(start:last));
    tempFFT(1) = [];
    n = length(tempFFT);
    tempFFT = abs(tempFFT(1:floor(n/2))).^2;
    height = length(tempFFT);
    myFFT = table('Size',[height,cz],'VariableTypes',varTypes,'VariableNames',varNames);
  
    %body forces and moments
    % get mean of when no jet is on
    ForceNoJet = mean(dataTable{start0:last0,[vars(1:6)]},1,"omitnan");
    
    BodyForce = dataTable{start:last,[vars(1:6)]};

    NetForce = BodyForce - ForceNoJet; % what will be written to Fx, Fy, ...
    
    TotalNetForce = sqrt(NetForce(:,1).^2+NetForce(:,2).^2+NetForce(:,3).^2);
    
    % compute FFTs
    NetForceFFT = fft(NetForce);
    NetForceFFT(1,:) =[]; % remove the first value (sum of the data)
    n = length(NetForceFFT);
    NetForceFFT_Power = abs(NetForceFFT(1:floor(n/2),:)).^2;
    
    % frequency range
    Freq = (1:n/2)/(n/2)*maxfreq;

    % total net force
    TotalNetForceFFT = fft(TotalNetForce);
    TotalNetForceFFT(1) = [];
    n = length(TotalNetForceFFT);
    TotalNetForceFFT_Power = abs(TotalNetForceFFT(1:floor(n/2))).^2;

    % Rotor Reaction Force
    FFTRotorThrust = fft(dataTable.RotorThrust(start:last));
    FFTRotorThrust(1) = [];
    n = length(FFTRotorThrust);
    RotorThrust_Power = abs(FFTRotorThrust(1:floor(n/2))).^2;

    % Rotor speed
    FFTRotorSpeed = fft(dataTable.("RotorSpeed-Frequency")(start:last));
    FFTRotorSpeed(1) =[];
    n = length(FFTRotorSpeed);
    RotorSpeed_Power = abs(FFTRotorSpeed(1:floor(n/2))).^2; 
    
    myFFT{:,varNames(1:6)} = NetForceFFT_Power;
    myFFT{:,"RotorThrust Power"} = RotorThrust_Power;
    myFFT{:,"RotorSpeed-Frequency Power"} = RotorSpeed_Power;
    myFFT{:,"Total Net Force Power"} = TotalNetForceFFT_Power;
    myFFT{:,"Frequency"} = Freq';
    myFFT{:,"WaterjetSpeed"} = testMatrix.(testMatrixVars{1})(testMatrixIndex)*ones(floor(n/2),1);
    myFFT{:,"TrimAngle"} = testMatrix.(testMatrixVars{2})(testMatrixIndex)*ones(floor(n/2),1);
  
    clear temp
    return
end