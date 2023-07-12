function waterjetForceFigureMaker(data)
%Inputs:
%data: data set containing all the necessary data for plotting
%index: Y-axis item
%label: Y-axis title
%barelabel: label w/o units
%length scale: used for non dimensionalization

%Ydata = (data, water depth, heading, steering, flow speed)

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

    % User selects which Item they would like on the y-axis
    yDataSelection = VariableNames;
    [indY,tf] = listdlg('ListString',string(VariableNames),'PromptString','Select y-Axis Variables');
    if ~tf
        fprintf("No selection made, exiting function");
        return;
    end
    % we don't want to plot Y = X so remove the IndY's from the indX
    % selection listdlg
    xVariableNames = VariableNames;
    xVariableNames(indY) = [];
    xDataSelection = xVariableNames;
    [indX,tf] = listdlg('ListString',string(xVariableNames),'SelectionMode','single','PromptString','Select x-Axis Variable');
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
    xData2Plot = xDataSelection(indX);

    % Modify yData to accept number of new elements for STDEV
    numberVariables = length(indY)+length(indX);
    indicies = [indY,indX];
    variables2Plot = [yData2Plot,xData2Plot];
    for i = 1:numberVariables
        yData(:,end+1) = {0}; % initialize
        newVarName = strcat("STDEV ",variables2Plot(i));
        VariableNames(end+1) = newVarName;
        stdevIndicies(i) = length(VariableNames);
    end
    yData.Properties.VariableNames = VariableNames;

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
    % Consolidate repeated trials, including TrimAngle, WaterjetSpeed, and
    % FlowSpeed
    if answerJSweep == "Bollard Pull"
        [GC,GR] = groupcounts(yData{:,"Waterjet Speed"});
        if any(GC>1)
            %there are repeats, identify the value
            indGC = find(GC>1);
            valueGR = GR(indGC); %this may be larger than 1 value
            %find corresponding row from yData3
            for a = 1:length(valueGR)
                rowsYData = find(yData{:,"Waterjet Speed"}==valueGR(a));
                temp = yData(rowsYData,:);
                yData(rowsYData,:) = []; %deletes the row
                mu = mean(temp{:,variables2Plot});
    
                if GC(indGC(a))>2
                    sigma = std(temp{:,variables2Plot});
                else
                    sigma = abs(temp{1,variables2Plot}-temp{2,variables2Plot});
                end
                newYDatarow = temp(1,:);
                for b = 1:length(indicies)
                    newYDatarow{1,indicies(b)} = mu(b); %should write to the correct place
                    newYDatarow{1,stdevIndicies(b)} = sigma(b);
                end
                yData(end+1,:) = newYDatarow;
            end
        end
        yData = sortrows(yData,"Waterjet Speed");
    else
        yData = sortrows(yData,xData2Plot);
    end
    % build figname string by concatinating ydata2plot strings
    for i = 1:length(indY)
        if i == 1
            figString = yData2Plot(i);
        else
            figString = strcat(figString," and ",yData2Plot(i));
        end
    end
    figname = strcat(figString," v. ", xData2Plot);
    % pull units from data for use in building labels
    units = data.WJ001{7:end,"Units"};
    units = [units(1:2); "RPM"; units(3:end)];
    units2plot = units(indY);
    figure("Name",figname,'units','normalized','OuterPosition',[0 0 1 1]); %makes full screen size
    for i = 1:length(indY)
        if i == 1
            ylabelString = strcat(yData2Plot(i)," (",units2plot(i),")");
        else
            ylabelString = strcat(ylabelString,"; ",yData2Plot(i)," (",units2plot(i),")");
        end
    end
    hold on
    for i = 1:length(indY)
        yneg = yData{:,stdevIndicies(i)};
        ypos = yneg;
        xneg = yData{:,stdevIndicies(end)}; % x err always last
        xpos = xneg;
        errorbar(yData{:,xData2Plot},yData{:,yData2Plot(i)},yneg,ypos,...
            xneg,xpos,"-s","MarkerFaceColor",myMonoChromeColorMap(i,:),...
            "MarkerEdgeColor",myColorMap(i,:),"Color",myMonoChromeColorMap(i,:),...
            "MarkerSize",10,'LineWidth',1,"DisplayName",yData2Plot(i));
    end
    ylabel(ylabelString);
    xUnits = units;
    xUnits(indY) = [];
    xlabel(strcat(xData2Plot," (",xUnits(indX),")"));
    title(figname);
    legend;
end

function myCloseReq(src,event)
    delete(src);
end