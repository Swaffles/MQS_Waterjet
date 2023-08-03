function waterjetFFTPlot(data)
%Inputs:
%data: data set containing all the necessary data for plotting

rho = 1000; %water density
gravity = 9.81; %acceleration due to gravity

%number of items to loop over
fields = fieldnames(data);
    for i=1:length(fields)
        yData(i,1) = data.(fields{i})(7,2);   % Rotor Thurst (N)
        yData(i,2) = data.(fields{i})(8,2);   % Rotor Speed (RPS)
        yData(i,3) = {yData{i,2}*60};         % Rotor Speed (RPM)
        yData(i,4) = data.(fields{i})(9,2);   % Water Depth (cm) 
        yData(i,5) = data.(fields{i})(10,2);  % Steering Angle
        temp1 = data.(fields{i}){11,2};
        yData(i,6) = {round((temp1/255)*100)}; % Motor Duty Cycle
        yData(i,7) = data.(fields{i})(12,2);  % Trim Angle
        temp2 = data.(fields{i}){13,2}/100;    % Flow Speed conversion
        yData(i,8) = {temp2};                  % Flow Speed (m/s)
        yData(i,9) = data.(fields{i})(14,2);  % Total Force (N)
        yData(i,10) = data.(fields{i})(15,2);  % Advance Coefficient J
        yData(i,11) = data.(fields{i})(16,2);  % Coefficient of thrust Kt
        yData(i,12) = data.(fields{i})(17,2); % Coefficient of thrust Kt 
        % rotor 
        %yData(i,13) = {0};                    % this will be where standard
        % deviation gets put in
    end
    VariableNames = [data.(fields{1}){7:17,1}]';
    VariableNames = [VariableNames(1:2) "Rotor Tachometer" VariableNames(3:end)];
    yData.Properties.VariableNames = VariableNames;
    clear temp1 temp2

    % UI Input
    % User selects Bollard pull or J Sweep
    answerJSweep = questdlg("What J condition?",...
        "J Condition","Bollard Pull", "J Sweep","Bollard Pull");
    switch answerJSweep
        case "Bollard Pull"
            % only the 0 J's remain
            indx = yData{:,"J"} == 0;
            yData = yData(indx,:);
        case "J Sweep"
            % 0 J's removed
            indx = yData{:,"J"} ~=0;
            yData = yData(indx,:);
        otherwise
            fprintf("Exiting function");
            return;
    end
    % User selects low submergence or normal
    answer = questdlg("What submergence condition","Submergence Condition",...
        "Normal","Low Submergence","Normal");
    switch answer
        case "Normal"
            % No zero trim angles
            indx = yData{:,"Trim Angle"} ~= 0;
            yData = yData(indx,:);
        case "Low Submergence"
            indx = yData{:,"Trim Angle"} == 0;
            yData = yData(indx,:);
        otherwise
           fprintf("Exiting function");
            return;
    end
    % User selects which Item(s) they would like to take FFT of
    yDataSelection = VariableNames;
    [indY,tf] = listdlg('ListString',string(VariableNames),'PromptString','Select y-Axis Variables');
    if ~tf
        fprintf("No selection made, exiting function");
        return;
    end
    % Graph can't handle more than four y-axis vairables. So we warn the
    % user and cut everything but the first four elements in indY. User
    % retains indX choice
    if length(indY)>4
        fig = uifigure;
        message = {'More than four (4) y-axis variables selected!','Cutting at element 4.'};
        uialert(fig,message,'Warning','Icon','warning',"CloseFcn",@myCloseReq);
        waitfor(fig)
        indY = indY(1:4);
    end
    yData2Plot = yDataSelection(indY);

    % Modify yData to accept number of new elements for STDEV
    numberVariables = length(indY);

    % edge colors
    myColorMap = [100/255 143/255 255/255;...
                  220/255 38/255 127/255;...
                  254/255 97/255 0/255;...
                  255/255 176/255 0/255]; %IBM color map  
    % face and line colors
    myMonoChromeColorMap = [136/255 163/255 230/255;...
                            194/255 72/255 131/255;...
                            229/255 116/255 46/255;...
                            230/255 172/255 46/255];
    % Outer loop to loop over each waterjet speed (motor counts)
    uniqueWaterJetSpeeds = unique(yData{:,"Waterjet Speed"});
    numWaterJetSpeeds = length(uniqueWaterJetSpeeds);
    for a = 1:numWaterJetSpeeds
        % fft of each selected variable
        inda = yData{:,"Waterjet Speed"} == uniqueWaterJetSpeeds(a);
        yData1 = yData(inda,:);
        for i = 1:numberVariables
            % build figname string by concatinating ydata2plot strings
            if i == 1
                figString = strcat("Power ",yData2Plot(i));
            else
                figString = strcat(figString," and Power ",yData2Plot(i));
            end
        end
        figname = strcat(figString," v. Frequency for motor duty cycle ",string(uniqueWaterJetSpeeds(a)),"%");
        figure("Name",figname,'units','normalized','OuterPosition',[0 0 1 1]); %makes full screen size
        hold on
        for k = 1:numberVariables
            y = fft(yData1{:,yData2Plot(k)});
            tempStr = strrep(yData2Plot(k),' ','_');
            y(1) = []; % remove the first value (sum of the data)
            n = length(y);
            yDataFFT.(tempStr).Power = abs(y(1:floor(n/2))).^2;
            maxfreq = 1000; % sampled at 2000
            yDataFFT.(tempStr).Freq = (1:n/2)/(n/2)*maxfreq;
            % plot it
            plot(yDataFFT.(tempStr).Freq,yDataFFT.(tempStr).Power,"Color",...
                myMonoChromeColorMap(k,:),'LineWidth',1,"DisplayName",...
                yData2Plot(k));
        end
        ylabel("Power");
        xlabel("Hz");
        title(figname);
        legend;
        hold off    
    end
    
    
    
    
end

function myCloseReq(src,event)
    delete(src);
end