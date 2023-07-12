function myTable = meanandstdev(dataTable,testMatrix,testMatrixVars,units,dataFileName,debug)
% This function takes in a table containting raw data from a labview
% testing log called dataTable, the corresponding testMatrix, the 
% dataFileNames, as well as the option to set a debugger. The function
% queries the user to select the bounds for valid data. The function then 
% outputs a table called myTable that contains the Quantity, the mean, 
% the STDDEV, and the units for each value from the input table and 
% specific values of the testMatrix.

% E.g.
% Input dataTable:
% Fx    Fy      Fz      ...
% 101   2.99
% 99    3.03
% .     .       .       ...
% .     .       .       ...

% Output myTable:
% Quantity      Mean        STDEV
% Fx            100         5
% Fy            3           0.04
% .             .           .
% .             .           .




     %locate the trial in the test matrix incase missaligned
    testMatrixIndex = find(strcmp(dataFileName,testMatrix.TrialName));
    if debug
        fprintf("Building myTable for %s\n", dataFileName);
        %fprintf("Using data from row %0.f of test matrix, corresponding to %s\n",...
        %        testMatrixIndex,testMatrix.TrialName(testMatrixIndex));
    end
    vars = fieldnames(dataTable);
    vars = vars(2:end-3); %drop last three table intrinsic properties
    vars(10) = []; %
    vars = [vars;testMatrixVars'];
    vars = strrep(vars,'_','-');
    sz = length(vars);
    varNames = ["Quantity","Mean","STDDEV","Units"];
    varTypes = ["string","double","double","string"];
    myTable = table('Size',[sz,4],'VariableTypes',varTypes,'VariableNames',varNames);
    myTable{:,end} = units;
  
    % have the user choose where steady state data is.
    for j=1:sz-length(testMatrixVars)
        plot(dataTable.(vars{j}));
    end
    [t,~] = ginput(2);
    start = round(t(1));
    last = round(t(2));
    %body forces and moments
    try
        meanBodyForce = mean(dataTable{start:last,[vars(1:6)]},1,'omitnan');

        stdBodyForce = std(dataTable{start:last,[vars(1:6)]},1,'omitnan');
    catch
        meanBodyForce = mean(dataTable{start:last,[vars(1:6)]},1,'omitnan');

        stdBodyForce = std(dataTable{start:last,[vars(1:6)]},1,'omitnan');
    end
    myTable{1:6,"Quantity"} = vars(1:6);
    myTable{1:6,["Mean","STDDEV"]} =...
           [meanBodyForce',stdBodyForce'];
    
    % water depth and port steering angle
    meanWaterDepth = mean(dataTable.WaterDepth,'omitnan');
    stdWaterDepth = std(dataTable.WaterDepth,'omitnan');
    % Rotor Reaction Force
    meanRotorThrust = mean(dataTable.RotorThrust,'omitnan');
    stdRotorThrust = std(dataTable.RotorThrust,'omitnan');
    % Steering Angle
    meanSteeringAngle = mean(dataTable.("SteeringAngle-PulseCount"),'omitnan');
    stdSteeringAngle = std(dataTable.("SteeringAngle-PulseCount"),'omitnan');
    % Rotor Speed
    meanRotorSpeed = mean(dataTable.("RotorSpeed-Frequency"),'omitnan');
    stdRotorSpeed = std(dataTable.("RotorSpeed-Frequency"),'omitnan');

    vars(7) = "Water Depth";
    vars(8) = "Rotor Thrust";
    vars(9) = "Steering Angle";
    vars(10) = "Rotor Speed";
    vars(11) = "Waterjet Speed";
    vars(12) = "Trim Angle";
    vars(13) = "Flow Speed";

    myTable{7:10,"Quantity"} = [vars(8);vars(10);vars(7);vars(9)];
    myTable{7:10,["Mean","STDDEV"]} =...
           [[meanRotorThrust;meanRotorSpeed;meanWaterDepth;meanSteeringAngle],...
           [stdRotorThrust;stdRotorSpeed;stdWaterDepth;stdSteeringAngle]];

    n = 11;
    for i=1:length(testMatrixVars)
        if i == 3
            %read both flow speeds and avg
            temp = testMatrix.Flow_meter_start(testMatrixIndex);
            temp(2) = testMatrix.Flow_meter_end(testMatrixIndex);
            %stddev by difference between start and end speed
            absDiff = abs(temp(1)-temp(2));
            myTable(n,:) = {vars{n},mean(temp),absDiff,units{n}};
        else
            myTable(n,:) = {vars{n},testMatrix.(testMatrixVars{i})(testMatrixIndex),...
                            0,units{n}};
        end
        n = n+1;
    end
%     try
%         myTable{end,4} = buoyantCorrection{end,2};
%     catch
%         fprintf("Error assigning H/D!\n");
%     end

    clear temp
    return
end