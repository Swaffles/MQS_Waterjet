% MQS Waterjet flume experiment data cruncher
% Objective: A quick data analysis tool for the waterjet tests.
% INPUTS: Test data (either TestLog## or WJ###) and the test matrix
% Operations: 1) Build a data table of values
% 2) Calculates the Mean and standard deviation for the data.
% 3) Create Figure plots of compiled data to form closed-water curves.

% note if you get a cell read error in meanandstdevARC it's because you have
% included a data item to read that has not been read from the excel sheet
% to the testMatrix table

clearvars -except homePath dataPath programPath
close all
clc

debug = true;

try 
    addpath(programPath)
catch
    programPath = uigetdir(pwd,'Select GitHub Folder');
    addpath(programPath);
end
try 
    addpath(dataPath)
catch
    % C:\Users\swafford\OneDrive - University of Iowa\MQS Waterjet Flume Test\Data
    dataPath = uigetdir(pwd,'Select University of Iowa\MQS\Waterjet Flume Test\Data');
end

try
    addpath(homePath)
    cd(homePath)
catch
    % C:\Users\swafford\OneDrive - University of Iowa\MQS Waterjet Flume Test
    homePath = uigetdir(pwd,'Select University of Iowa \MQS Waterjet Flume Test');
    cd(homePath)
end

%savePath = uigetdir(homePath,'Select Folder for Image file Saving');

if ~isfile("WaterjetData.mat") || debug
    %if running for first time or if debug is on
    noData = true;
    existingData = [];
    try
        load("waterjetTimeData.mat");
        TimeData = Data;
    catch
        fprintf("No time data :(\n");
    end
else
    load("FFTwaterjetData.mat");
    FFTData = Data;
    load("waterjetTimeData.mat");
    TimeData = Data;
    load('WaterjetData.mat');
    fprintf("Loading saved data from disc...\n");
    existingData = fieldnames(Data);
    noData = false; 
end

%gather all the data
cd (dataPath);

%read the exel first
cd ..
testMatrixFile = dir('Waterjet Test Matrix.xlsx');
% Specify sheet and range
opts = spreadsheetImportOptions("NumVariables", 16);

opts.Sheet = "Test Matrix";
opts.DataRange = "A4:P95";

% Specify column names and types
opts.VariableNames = ["TrialName", "WaterjetSpeed", "J", "WaterSpeed",...
    "WaterDepth", "TrimAngle", "CGx", "Position",...
    "PumpDutyCycle", "WaterjetAngle", "Duration", "Flow_meter_start",...
    "Flow_meter_end","Ultrasonic", "GoPro", "Comments"];
opts.VariableTypes = ["string", "double", "double", "double",...
    "double", "double", "double", "double",...
    "double", "double", "double", "double",...
    "double", "double", "double", "string"];

% Specify variable properties
opts = setvaropts(opts, "TrialName", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["TrialName", "Comments"], "EmptyFieldRule", "auto");
fprintf('Reading in test matrix file...');
testMatrix = readtable(testMatrixFile.name,opts);

cd (dataPath);
dataFiles = dir('WJ*');
dataFileNames = {dataFiles.name};
%read dataFiles into Table for processing
testMatrixVars = ["WaterjetSpeed","TrimAngle","FlowSpeed"];
units=["N","N","N","N-m","N-m","N-m","N","rps","cm","deg","% duty cycle","deg","cm/s"]';
if length(dataFiles)>length(existingData)
    saveMe = true;
    for ind = 1:length(dataFiles)
        if noData
           [~,~,fExt] = fileparts(dataFileNames{ind});
           if strcmp('.mat',fExt)
               fprintf("Building data table entry for %s\n", dataFileNames{ind});
               T = load(dataFileNames{ind});
               fieldName = fieldnames(T);
               T = T.(fieldName{1});
               temp = T.Properties.VariableNames;
               temp2 = strrep(temp,'_','-');
               T = renamevars(T,temp,temp2);
               f1 = figure('Name',strcat(dataFileNames{ind},' T Selection'),'units','normalized','OuterPosition',[0 0 1 1]);
               hold on
               % have the user choose where steady state data is.
               vars = fieldnames(T);
               vars = vars(2:end-3); %drop last three table intrinsic properties
               vars(10) = []; %
               sz = length(vars);
               try
                   t = TimeData.(fieldName{1});
               catch
                   for j=1:sz
                       plot(T.(vars{j}));
                   end
                   [t,~] = ginput(4);
               end
               start0 = round(t(1)); % user chooses zero data start
               last0 = round(t(2));  % user chooses zero data end
               start = round(t(3));  % user chooses steady state start
               last = round(t(4));   % user chooses steady state end 
               % if user picks a point outside of the table, usually the last data
               % point, select the last data point as the last point
               if last > size(T,1)
                   last = size(T,1); % +1 because we deleted the first element in tempFFT
               end
               time = [start0,last0,start,last];
               Data.(fieldName{1}) = meanandstdev(T,testMatrix,testMatrixVars,units,fieldName{1},time,debug);
               FFTData.(fieldName{1}) = fftWaterJet(T,testMatrix,testMatrixVars,fieldName{1},time,debug);
               PSDdata.(fieldName{1}) = psdWaterJet(T,testMatrix,testMatrixVars,fieldName{1},time,debug);
               Time.(fieldName{1}) = time;
               hold off
               close all
           else
               opts = detectImportOptions(dataFileNames{ind});
               opts.VariableNamingRule = 'modify';
               opts.VariableNamesLine = 7;
               opts.Delimiter = '\t';
               opts.VariableUnitsLine = 5;
               opts.DataLines = 8;
               % Specify file level properties
               opts.ExtraColumnsRule = "ignore";
               opts.EmptyLineRule = "read";
               fprintf("Building data table entry for %s\n", dataFileNames{ind});
               T = readtable(dataFileNames{ind},opts);
               temp = T.Properties.VariableNames;
               temp2 = strrep(temp,'_','-');
               T = renamevars(T,temp,temp2);
               f1 = figure('Name',strcat(dataFileNames{ind},' T Selection'),'units','normalized','OuterPosition',[0 0 1 1]);
               hold on
               % have the user choose where steady state data is.
               vars = fieldnames(T);
               vars = vars(2:end-3); %drop last three table intrinsic properties
               vars(10) = []; %
               sz = length(vars);
               try
                   t = TimeData.(dataFileNames{ind});
               catch
                   for j=1:sz
                       plot(T.(vars{j}));
                   end
                   [t,~] = ginput(4);
               end
               start0 = round(t(1)); % user chooses zero data start
               last0 = round(t(2));  % user chooses zero data end
               start = round(t(3));  % user chooses steady state start
               last = round(t(4));   % user chooses steady state end
               % if user picks a point outside of the table, usually the last data
               % point, select the last data point as the last point
               if last > size(T,1)
                   last = size(T,1); % +1 because we deleted the first element in tempFFT
               end
               time = [start0,last0,start,last];
               Data.(dataFileNames{ind}) = meanandstdev(T,testMatrix,testMatrixVars,units,dataFileNames{ind},time,debug);
               FFTData.(dataFileNames{ind}) = fftWaterJet(T,testMatrix,testMatrixVars,dataFileNames{ind},time,debug);
               PSDdata.(dataFileNames{ind}) = psdWaterJet(T,testMatrix,testMatrixVars,dataFileNames{ind},time,debug);
               Time.(dataFileNames{ind}) = time;
               hold off
               close all
           end
           saveData(Data,"WaterjetData.mat",homePath,dataPath);
           saveData(FFTData,"FFTwaterjetData.mat",homePath,dataPath);
           saveData(Time,"waterjetTimeData.mat",homePath,dataPath);
        else
            %if the file isn't already a part of T add it
            if ind>length(existingData)
               [~,~,fExt] = fileparts(dataFileNames{ind});
               if strcmp('.mat',fExt)
                   fprintf("Building data table entry for %s\n", dataFileNames{ind});
                   T = load(dataFileNames{ind});
                   fieldName = fieldnames(T);
                   T = T.(fieldName{1});
                   temp = T.Properties.VariableNames;
                   temp2 = strrep(temp,'_','-');
                   T = renamevars(T,temp,temp2);
                   f2 = figure('Name',strcat(dataFileNames{ind},' T Selection'),'units','normalized','OuterPosition',[0 0 1 1]);
                   hold on
                   % have the user choose where steady state data is.
                   vars = fieldnames(T);
                   vars = vars(2:end-3); %drop last three table intrinsic properties
                   vars(10) = []; %
                   sz = length(vars);
                   try
                       t = TimeData.(fieldName{1});
                   catch
                       for j=1:sz
                           plot(T.(vars{j}));
                       end
                       [t,~] = ginput(4);
                   end
                   start0 = round(t(1)); % user chooses zero data start
                   last0 = round(t(2));  % user chooses zero data end
                   start = round(t(3));  % user chooses steady state start
                   last = round(t(4));   % user chooses steady state end
                   % if user picks a point outside of the table, usually the last data
                   % point, select the last data point as the last point
                   if last > size(T,1)
                       last = size(T,1); % +1 because we deleted the first element in tempFFT
                   end
                   time = [start0,last0,start,last];
                   Data.(fieldName{1}) = meanandstdev(T,testMatrix,testMatrixVars,units,fieldName{1},time,debug);
                   FFTData.(fieldName{1}) = fftWaterJet(T,testMatrix,testMatrixVars,fieldName{1},time,debug);
                   PSDdata.(fieldName{1}) = psdWaterJet(T,testMatrix,testMatrixVars,fieldName{1},time,debug);
                   Time.(fieldName{1}) = time;
                   hold off
                   close all
               else
                   opts = detectImportOptions(dataFileNames{ind});
                   opts.VariableNamingRule = 'modify';
                   opts.VariableNamesLine = 7;
                   opts.Delimiter = '\t';
                   opts.VariableUnitsLine = 5;
                   opts.DataLines = 8;
                   % Specify file level properties
                   opts.ExtraColumnsRule = "ignore";
                   opts.EmptyLineRule = "read";
                   fprintf("Building data table entry for %s\n", dataFileNames{ind});
                   T = readtable(dataFileNames{ind},opts);
                   temp = T.Properties.VariableNames;
                   temp2 = strrep(temp,'_','-');
                   T = renamevars(T,temp,temp2);
                   f2 = figure('Name',strcat(dataFileNames{ind},' T Selection'),'units','normalized','OuterPosition',[0 0 1 1]);
                   hold on
                   % have the user choose where steady state data is.
                   vars = fieldnames(T);
                   vars = vars(2:end-3); %drop last three table intrinsic properties
                   vars(10) = []; %
                   sz = length(vars);
                   try
                       t = TimeData.(dataFileNames{ind});
                   catch
                       for j=1:sz
                           plot(T.(vars{j}));
                       end
                       [t,~] = ginput(4);
                   end
                   start0 = round(t(1)); % user chooses zero data start
                   last0 = round(t(2));  % user chooses zero data end
                   start = round(t(3));  % user chooses steady state start
                   last = round(t(4));   % user chooses steady state end
                   % if user picks a point outside of the table, usually the last data
                   % point, select the last data point as the last point
                   if last > size(T,1)
                       last = size(T,1); % +1 because we deleted the first element in tempFFT
                   end
                   time = [start0,last0,start,last];
                   Data.(dataFileNames{ind}) = meanandstdev(T,testMatrix,testMatrixVars,units,dataFileNames{ind},time,debug);
                   FFTData.(dataFileNames{ind}) = fftWaterJet(T,testMatrix,testMatrixVars,dataFileNames{ind},time,debug);
                   PSDdata.(dataFileNames{ind}) = psdWaterJet(T,testMatrix,testMatrixVars,dataFileNames{ind},time,debug);
                   Time.(dataFileNames{ind}) = time;
                   hold off
                   close all
               end
            saveData(Data,"WaterjetData.mat",homePath,dataPath);
            saveData(FFTData,"FFTwaterjetData.mat",homePath,dataPath);
            saveData(Time,"waterjetTimeData.mat",homePath,dataPath);
            end
        end
        clear T
    end
else
    %diff = length(dataFiles) - length(existingData)
    fprintf("No new files\n");
    saveMe = false;
end


%%
%Individual file key metrics

% Calculate Total Force
fields = fieldnames(Data);
for ind = 1:length(fields)
    Fx = Data.(fields{ind}){1,2};
    Fy = Data.(fields{ind}){2,2};
    Fz = Data.(fields{ind}){3,2};
    totalForce = sqrt(Fx^2+Fy^2+Fz^2);
    % we need to determine whether this force is thrust or drag. If drag
    % do not calculate a thrust coefficient. Fz is positive downwards!
    if Fx<0
        totalForce = -1*totalForce;
    end

    sFx = Data.(fields{ind}){1,3};
    sFy = Data.(fields{ind}){2,3};
    sFz = Data.(fields{ind}){3,3};
    sTotalForce = sqrt(sFx^2+sFy^2+sFz^2);
    Data.(fields{ind}){14,:} = ["Total Force",totalForce,sTotalForce,"N"];
end
% Calculate Advance Coefficient J
rotorDiameter = 50.8/1000; % 50.8 mm rotor
for ind = 1:length(fields)
    U = Data.(fields{ind}){13,2}; % flow speed in cm/s
    U = U/100; % convert to meters
    N = Data.(fields{ind}){8,2}; % rotor speed in rps
    J = U/(N*rotorDiameter);

    sU = Data.(fields{ind}){13,3};
    sU = sU/100;
    sN = Data.(fields{ind}){8,3};
    thetaU = 1/(N*rotorDiameter);
    thetaN = -U/(N^2*rotorDiameter);
    sJ = sqrt(thetaU^2*sU^2+thetaN^2*sN^2);

    Data.(fields{ind}){15,:} = ["J",J,sJ,"--"];
end
% Calculate KT total
rho = 1000;
for ind = 1:length(fields)
    Tm = Data.(fields{ind}){14,2};

        N = Data.(fields{ind}){8,2};
        KT = Tm/(rho*rotorDiameter^4*N^2);
    
        sTm = Data.(fields{ind}){14,3};
        sN = Data.(fields{ind}){8,3};
        thetaTm = 1/(rho*rotorDiameter^4*N^2);
        thetaN = (-2*Tm)/(rho*rotorDiameter^4*N^3);
        sKT = sqrt(thetaTm^2*sTm^2+thetaN^2*sN^2);
    
        Data.(fields{ind}){16,:} = ["Kt",KT,sKT,"--"];
end
% Calculate KT rotor
for ind = 1:length(fields)
    RTm = Data.(fields{ind}){7,2};
    N = Data.(fields{ind}){8,2};
    KTr = RTm/(rho*rotorDiameter^4*N^2);

    sRTm = Data.(fields{ind}){7,3};
    sN = Data.(fields{ind}){8,3};
    thetaRTm = 1/(rho*rotorDiameter^4*N^2);
    thetaN = (-2*RTm)/(rho*rotorDiameter^4*N^3);
    sKTr = sqrt(thetaRTm^2*sRTm^2+thetaN^2*sN^2);

    Data.(fields{ind}){17,:} = ["Ktr",KTr,sKTr,"--"];
end
FFTfields = fieldnames(FFTData);
for ind = 1:length(FFTfields)
    FFTData.(FFTfields{ind}){:,"J"} = Data.(fields{ind}).Mean(15);
end
PSDfields = fieldnames(PSDdata);
for ind = 1:length(PSDfields)
    PSDdata.(PSDfields{ind}){:,"J"} = Data.(fields{ind}).Mean(15);
end
saveData(Data,"WaterjetData.mat",homePath,dataPath);
saveData(FFTData,"FFTwaterjetData.mat",homePath,dataPath);
saveData(PSDdata,"PSDwaterjetData.mat",homePath,dataPath);
clearvars -except Data FFTData PSDdata Time vars homePath dataPath programPath savePath debug 
%% Graph Total Force and Rotor Force v. Shaft RPM Bollard Pull
if debug
    close all
end
cd(homePath);
waterjetForceFigureMaker(Data);
hold off

%% FFT
if debug
    close all
end
cd(homePath);
waterjetFFTPlot(FFTData);
hold off

%% PSD
if debug
    close all
end
cd(homePath)
waterjetPSDPlot(PSDdata);
hold off
%% SAVE Function
function saveData(Data,filenameString,homePath,dataPath)
    cd (homePath);
    %save in homePath
    text = strcat("Saving tables to ",filenameString);
    fprintf(text);
    fprintf(" in %s\n", homePath);
    save(filenameString,"Data",'-v7.3');
    cd (dataPath);
end