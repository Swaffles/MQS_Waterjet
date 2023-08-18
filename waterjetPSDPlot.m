function waterjetPSDPlot(data)
%Inputs:
%data: data set containing all the necessary data for plotting

rho = 1000; %water density
gravity = 9.81; %acceleration due to gravity

%number of items to loop over
fields = fieldnames(data);
tableSize = [length(fields),width(data.WJ001)+1];
tableVars = ["Trial",data.WJ001.Properties.VariableNames];
tableVarTypes = ["string",varfun(@class,data.WJ001,"OutputFormat","cell")];
yData = table('Size',tableSize,'VariableTypes',tableVarTypes,'VariableNames',tableVars);
    for i=1:length(fields)
        yData(i,1) = fields(i);                     % trial name
        yData(i,2:end) = head(data.(fields{i}),1);  % 1st line of data
    end
    % VariableNames = [data.(fields{1}){7:17,1}]';
    % VariableNames = [VariableNames(1:2) "Rotor Tachometer" VariableNames(3:end)];

    VariableNames = yData.Properties.VariableNames;

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
            indx = yData{:,"TrimAngle"} ~= 0;
            yData = yData(indx,:);
        case "Low Submergence"
            indx = yData{:,"TrimAngle"} == 0;
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
    % myColorMap = [100/255 143/255 255/255;...
    %               220/255 38/255 127/255;...
    %               254/255 97/255 0/255;...
    %               255/255 176/255 0/255]; %IBM color map  
    % face and line colors
    myMonoChromeColorMap = [136/255 163/255 230/255;... % Light blue
                            120/255 94/255 240/255;...  % Purple
                            194/255 72/255 131/255;...  % Pink
                            229/255 116/255 46/255;...  % Orange
                            230/255 172/255 46/255];    % Yellow
    myLineStyles = ["-","--",":","-."];
    % Outer loop to loop over each waterjet speed (motor counts)
    uniqueWaterJetSpeeds = unique(yData{:,"WaterjetSpeed"});
    uniqueWaterJetSpeeds = sort(uniqueWaterJetSpeeds,'descend');
    numWaterJetSpeeds = length(uniqueWaterJetSpeeds);
    for i = 1:numberVariables
        % build figname string by concatinating ydata2plot strings
        if i == 1
            figString = yData2Plot(i);
        else
            figString = strcat(figString," and ",yData2Plot(i));
        end
    end
    figname = strcat(figString," v. Frequency");
    figure("Name",figname,'units','normalized','OuterPosition',[0 0 1 1]); %makes full screen size
    tiledlayout(numWaterJetSpeeds,1);
    for a = 1:numWaterJetSpeeds
        % fft of each selected variable
        inda = yData{:,"WaterjetSpeed"} == uniqueWaterJetSpeeds(a);
        items = fields(inda);
        nexttile;
        tileTitle = strcat("Motor Duty Cycle ",...
        string(round((uniqueWaterJetSpeeds(a)/255)*100)),"%");
        hold on
        for k = 1:numberVariables
            % select the data to plot
            colorCount = 1;
            styleCount = 1;
            for j = 1:length(items)
                % plot it
                Freq = data.(items{j}).Frequency;
                Power = data.(items{j}).(yData2Plot{k});
                semilogy(Freq,Power,"Color",myMonoChromeColorMap(colorCount,:),...
                    'LineWidth',1,"LineStyle",myLineStyles(styleCount),"DisplayName",strcat(items(j)," ",yData2Plot(k)));
                grid on
                [~,locs] = findpeaks(Power,"MinPeakHeight",1.5e7,"NPeaks",6,"Threshold",1.5e7,"MinPeakProminence",1e9*0.02); % limit to 6 peaks
                sprintf('Peak locations for %s in %s:',items{j},yData2Plot{k});
                disp(Freq(locs));
                % need to check if we are out of colors
                if mod(j,5) == 0
                    colorCount = 1;
                    styleCount = styleCount+1;
                    if styleCount == 5
                        fprintf("Unique line styles and colors exceeded!\n");
                    end
                else
                    colorCount = colorCount+1;
                end
            end
        end
        ylabel("Power");
        xlabel("Hz");
        title(tileTitle);
        legend('Location','bestoutside');
        hold off
        grid off
    end
    
    
    
    
end

function myCloseReq(src,event)
    delete(src);
end