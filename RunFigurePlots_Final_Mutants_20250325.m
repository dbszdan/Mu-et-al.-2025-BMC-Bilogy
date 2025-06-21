baseStartingVolume = 9.6206E-17;

initialMolarityDelta = 1.0;
exoTimeDelay = 5.0;
initialk5Max = [1.5];
initialfEis = [0.05];
initial_eps_Ca = [0.01];
initialCaOut = 2.5E-3*1000;
colorList = ["blue" "red"];

gap_height = 4.5e-6*1E6;
small_r = gap_height*0.5;


% reference case
%[maxEpsValue, plotArray, t] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1), baseStartingVolume, initialfEis(1), initial_eps_Ca(1), initialCaOut, true, true, false);
[maxEpsValue, plotArray, t]  = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1), baseStartingVolume, initialfEis(1), initial_eps_Ca(1), initialCaOut, true, false, false);
%pause(1000);

% pil1D
[maxEpsValue_pil1D, plotArray_pil1D, t_pil1D] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1), baseStartingVolume, 0.0, initial_eps_Ca(1), initialCaOut, true, true, false);

% tricalbin delta
[maxEpsValue_tricalbD, plotArray_tricalbD, t_tricalbD] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, 0.0, baseStartingVolume, initialfEis(1), initial_eps_Ca(1), initialCaOut, true, true, false);

% no exocytosis, BFA
[maxEpsValue_exoD, plotArray_exoD, t_exoD] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1), baseStartingVolume, initialfEis(1), initial_eps_Ca(1), initialCaOut, false, true, false);
%pause(1000);

% no external Ca
[maxEpsValue_extCaD, plotArray_extCaD, t_extCaD] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1), baseStartingVolume, initialfEis(1), initial_eps_Ca(1), 0.0, true, true, false);

% tricalbin overexpression
%2.0x
[maxEpsValue_tricalbO, plotArray_tricalbO, t_tricalbO] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1)*2.0, baseStartingVolume, initialfEis(1), initial_eps_Ca(1), initialCaOut, true, true, false);

% no endocytosis, LatA
[maxEpsValue_endoD, plotArray_endoD, t_endoD] = AreaBalancePassValues_Final_20250325(initialMolarityDelta, exoTimeDelay, initialk5Max(1), baseStartingVolume, initialfEis(1), initial_eps_Ca(1), initialCaOut, true, false, false);
%pause(1000)

epsThres = 0.06;
epsThres = 0.056;
indexBurst = length(plotArray(:,12));
foundIndexBurst = false;

indexBurstMutants = [length(plotArray(:,12)) length(plotArray(:,12)) length(plotArray(:,12)) length(plotArray(:,12)) length(plotArray(:,12))  length(plotArray(:,12))];
foundIndexBurstMutants = [false false false false false false];

indexThreshold = [-1];
indexThreshold_pil1D = [-1];
indexThreshold_tricalbD = [-1];
indexThreshold_exoD = [-1];
indexThreshold_extCaD = [-1];
indexThreshold_tricalbO = [-1];
indexThreshold_endoD = [-1];

AreaByExo = zeros(1, length(t));
AreaByEndo = zeros(1, length(t));
netAreaByEisosome = zeros(1, length(t));
AreaByTricalibin = zeros(1, length(t));


for i = 1:length(plotArray(:,12))
    curr_eps_value = plotArray(i,12);

    if curr_eps_value > epsThres
        if foundIndexBurst == false
            indexBurst = i;
            foundIndexBurst = true;
        end

        

        if indexThreshold(1,1) == -1
            indexThreshold(1,1) = i;
        elseif mod(length(indexThreshold(1,:)), 2) == 0
            indexThreshold(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold(1,:)), 2) ~= 0 && indexThreshold(1,1) ~= -1
            indexThreshold(1,end+1) = i;
        end
    end
    
    %bad way to do this...
    if plotArray_pil1D(i,12) > epsThres 
        if foundIndexBurstMutants(1) == false
            indexBurstMutants(1) = i;
            foundIndexBurstMutants(1) = true;
        end

        if indexThreshold_pil1D(1,1) == -1
            indexThreshold_pil1D(1,1) = i;
        elseif mod(length(indexThreshold_pil1D(1,:)), 2) == 0
            indexThreshold_pil1D(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold_pil1D(1,:)), 2) ~= 0 && indexThreshold_pil1D(1,1) ~= -1
            indexThreshold_pil1D(1,end+1) = i;
        end
    end

    if plotArray_tricalbD(i,12) > epsThres 
        if foundIndexBurstMutants(2) == false
            indexBurstMutants(2) = i;
            foundIndexBurstMutants(2) = true;
        end

        if indexThreshold_tricalbD(1,1) == -1
            indexThreshold_tricalbD(1,1) = i;
        elseif mod(length(indexThreshold_tricalbD(1,:)), 2) == 0
            indexThreshold_tricalbD(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold_tricalbD(1,:)), 2) ~= 0 && indexThreshold_tricalbD(1,1) ~= -1
            indexThreshold_tricalbD(1,end+1) = i;
        end
    end


    if plotArray_exoD(i,12) > epsThres 
        if foundIndexBurstMutants(3) == false
            indexBurstMutants(3) = i;
            foundIndexBurstMutants(3) = true;
        end

        if indexThreshold_exoD(1,1) == -1
            indexThreshold_exoD(1,1) = i;
        elseif mod(length(indexThreshold_exoD(1,:)), 2) == 0
            indexThreshold_exoD(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold_exoD(1,:)), 2) ~= 0 && indexThreshold_exoD(1,1) ~= -1
            indexThreshold_exoD(1,end+1) = i;
        end
    end

    if plotArray_extCaD(i,12) > epsThres 
        if foundIndexBurstMutants(4) == false
            indexBurstMutants(4) = i;
            foundIndexBurstMutants(4) = true;
        end

        if indexThreshold_extCaD(1,1) == -1
            indexThreshold_extCaD(1,1) = i;
        elseif mod(length(indexThreshold_extCaD(1,:)), 2) == 0
            indexThreshold_extCaD(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold_extCaD(1,:)), 2) ~= 0 && indexThreshold_extCaD(1,1) ~= -1
            indexThreshold_extCaD(1,end+1) = i;
        end
    end

    if plotArray_tricalbO(i,12) > epsThres 
        if foundIndexBurstMutants(5) == false
            indexBurstMutants(5) = i;
            foundIndexBurstMutants(5) = true;
        end

        if indexThreshold_tricalbO(1,1) == -1
            indexThreshold_tricalbO(1,1) = i;
        elseif mod(length(indexThreshold_tricalbO(1,:)), 2) == 0
            indexThreshold_tricalbO(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold_tricalbO(1,:)), 2) ~= 0 && indexThreshold_tricalbO(1,1) ~= -1
            indexThreshold_tricalbO(1,end+1) = i;
        end
    end
    
    if plotArray_endoD(i,12) > epsThres 
        if foundIndexBurstMutants(6) == false
            indexBurstMutants(6) = i;
            foundIndexBurstMutants(6) = true;
        end

        if indexThreshold_endoD(1,1) == -1
            indexThreshold_endoD(1,1) = i;
        elseif mod(length(indexThreshold_endoD(1,:)), 2) == 0
            indexThreshold_endoD(1,end+1) = i; 
        end
    else
         if mod(length(indexThreshold_endoD(1,:)), 2) ~= 0 && indexThreshold_endoD(1,1) ~= -1
            indexThreshold_endoD(1,end+1) = i;
        end
    end
    
    AreaByExo(i) = plotArray(i,8)/60.0*(t(2)-t(1));
    AreaByEndo(i) = - plotArray(i,9)/60.0*(t(2)-t(1));
    netAreaByEisosome(i) = plotArray(i,10)/60.0*(t(2)-t(1)) - plotArray(i,11)/60.0*(t(2)-t(1));
    AreaByTricalibin(i) = plotArray(i,21)/60.0*(t(2)-t(1));
    
    if i > 1
        AreaByExo(i) = AreaByExo(i) + AreaByExo(i-1);
        AreaByEndo(i) = AreaByEndo(i) + AreaByEndo(i-1);
        netAreaByEisosome(i) = netAreaByEisosome(i) + netAreaByEisosome(i-1);
        AreaByTricalibin(i) = AreaByTricalibin(i) + AreaByTricalibin(i-1);
    end
    
end


%disp(indexThreshold_extCaD)

Area_values = 2.0*pi^2*small_r.*plotArray(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray(:,1) .* plotArray(:,1);

Area_values_normalized = Area_values/Area_values(1);



figure1 = figure(1);
clf()
%hold on

axes1 = axes('Parent',figure1,...
    'Position',[0.143214281988995 0.154603170558574 0.732500006865177 0.685873019917617]);
hold(axes1,'on');


colormap('jet')

yyaxis left
scatter(t(1:end), Area_values_normalized(1:end), 30, plotArray(:,12),"filled")
ylim([1.0, 1.3]);

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.5 0.5 0.5];

ax.FontSize = 14;

cb = colorbar('northoutside');
%a = get(cb);
%a = a.Position;
set(cb,'Position', [0.155357142857142 0.866666666666667 0.380357142857144 0.0333333324818564], 'FontSize',10)

%set(axes1,'CLim',[0 0.1],'FontSize',14);
%colorbar('northoutside', axes1,'Position',...
%    [0.155357142857142 0.866666666666667 0.380357142857144 0.0333333324818564],...
%    'FontSize',10);

clim([0.0, 0.10]);
ylabel(cb, 'Strain', 'FontSize', 10)

%{
max_eps = max(plotArray(:,12));
min_eps = min(plotArray(:,12));

cmap_fidelity = 1000;
curr_cmap = jet(cmap_fidelity);


for i = 1:length(t)-1
    % Plot each segment with a color based on the third array
    curr_index = (plotArray(i,12)-min_eps)/max_eps*cmap_fidelity-1)+1;
    disp(i)
    disp(min_eps)
    disp(plotArray(i,12))
    disp(curr_index)
    curr_color = curr_cmap(curr_index,:);
    %disp(curr_color)
    line([t(i), t(i+1)], [Area_values_normalized(i), Area_values_normalized(i+1)], 'Color', curr_color, 'LineWidth', 2);
    hold on;
end
%p = plot(t(1:end), Area_values_normalized(1:end), 'Color', cd);
%}


%plot(t(indexBurst:end), Area_values_normalized(indexBurst:end), 'b--');
ylabel('Relative surface area', 'FontSize', 20)

yyaxis right
plot(t, plotArray(:,20), 'Color', [0.5 0.5 0.5]);
ylabel('Ca level', 'FontSize', 20)

ylim([0.0, 2.5]);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

yyaxis left

if indexThreshold(1) ~= -1
    for i = 1:1:length(indexThreshold)
        plot([80-40,150-40], [Area_values_normalized(indexThreshold(i)), Area_values_normalized(indexThreshold(i))], 'Color', 'black', 'LineStyle','--', 'LineWidth',1.5)
    end
end

xlabel('Time (s)', 'FontSize', 20)


box on


figure(2)
clf()
hold on
plot(t, plotArray(:, 12), 'b');
plot([t(1) t(end)], [epsThres epsThres], 'r--');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

%disp(plotArray(:, 12))

ylim([-1E-5 0.075])

xlabel('Time (s)', 'FontSize', 20)
ylabel('Strain', 'FontSize', 20)

ax = gca();
ax.FontSize = 14;
box on

figure(3)
clf()
hold on

ax = gca();
ax.FontSize = 14;

plot([t(1) t(end)], [0.0, 0.0], '--', 'Color', [.7 .7 .7])
plot(t, AreaByExo, 'b');
plot(t, AreaByEndo, 'r')
plot(t, netAreaByEisosome, 'g');
plot(t, AreaByTricalibin, 'Color', [0.4940 0.1840 0.5560]);
%plot(t, AreaByExo+AreaByEndo+netAreaByEisosome+AreaByTricalibin, 'k');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);


xlabel('Time (s)', 'FontSize', 20)
ylabel('Net surface area addition (\mum^2)', 'FontSize', 20)

box on

ylim([-10 35])



figure(4)
clf()
hold on
plot(t(1:indexBurst), Area_values_normalized(1:indexBurst), 'k', 'LineWidth', 2);
plot(t(indexBurst:end), Area_values_normalized(indexBurst:end), 'k--', 'LineWidth', 2);
plot(t(indexBurst), Area_values_normalized(indexBurst), 'o', 'Color', 'k', 'LineWidth', 2);
plot([t(1) t(end)], [Area_values_normalized(indexBurst) Area_values_normalized(indexBurst)], 'Color', 'k', 'LineWidth', 1);

ax = gca();
ax.FontSize = 14;

mycolors = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 200 200 200];
mycolors = mycolors./255.0;

%pil1D
Area_values_mutant = 2.0*pi^2*small_r.*plotArray_pil1D(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray_pil1D(:,1) .* plotArray_pil1D(:,1);
Area_values_normalized_pil1D = Area_values_mutant ./ Area_values_mutant(1);

plot(t(1:indexBurstMutants(1)), Area_values_normalized_pil1D(1:indexBurstMutants(1)), 'Color', mycolors(1,:), 'LineWidth', 2);
plot(t(indexBurstMutants(1):end), Area_values_normalized_pil1D(indexBurstMutants(1):end), '--', 'Color', mycolors(1,:), 'LineWidth', 2);
plot(t(indexBurstMutants(1)), Area_values_normalized_pil1D(indexBurstMutants(1)), 'o', 'Color', mycolors(1,:), 'LineWidth', 2);
plot([t(1) t(end)], [Area_values_normalized_pil1D(indexBurstMutants(1)) Area_values_normalized_pil1D(indexBurstMutants(1))], 'Color', mycolors(1,:), 'LineWidth', 1);



%tricalbin delta
Area_values_mutant = 2.0*pi^2*small_r.*plotArray_tricalbD(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray_tricalbD(:,1) .* plotArray_tricalbD(:,1);
Area_values_normalized_tricalbD = Area_values_mutant ./ Area_values_mutant(1);

plot(t(1:indexBurstMutants(2)), Area_values_normalized_tricalbD(1:indexBurstMutants(2)), 'Color', mycolors(2,:), 'LineWidth', 2);
plot(t(indexBurstMutants(2):end), Area_values_normalized_tricalbD(indexBurstMutants(2):end), '--', 'Color', mycolors(2,:), 'LineWidth', 2);
plot(t(indexBurstMutants(2)), Area_values_normalized_tricalbD(indexBurstMutants(2)), 'o', 'Color', mycolors(2,:), 'LineWidth', 2);
plot([t(1) t(end)], [Area_values_normalized_tricalbD(indexBurstMutants(2)) Area_values_normalized_tricalbD(indexBurstMutants(2))], 'Color', mycolors(2,:), 'LineWidth', 1);

%no exocytosis
Area_values_mutant = 2.0*pi^2*small_r.*plotArray_exoD(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray_exoD(:,1) .* plotArray_exoD(:,1);
Area_values_normalized_exoD = Area_values_mutant ./ Area_values_mutant(1);

plot(t(1:indexBurstMutants(3)), Area_values_normalized_exoD(1:indexBurstMutants(3)), 'Color', mycolors(3,:), 'LineWidth', 2);
plot(t(indexBurstMutants(3):end), Area_values_normalized_exoD(indexBurstMutants(3):end), '--', 'Color', mycolors(3,:), 'LineWidth', 2);
plot(t(indexBurstMutants(3)), Area_values_normalized_exoD(indexBurstMutants(3)), 'o', 'Color', mycolors(3,:), 'LineWidth', 2);
plot([t(1) t(end)], [Area_values_normalized_exoD(indexBurstMutants(3)) Area_values_normalized_exoD(indexBurstMutants(3))], 'Color', mycolors(3,:), 'LineWidth', 1);

% no external Ca
Area_values_mutant = 2.0*pi^2*small_r.*plotArray_extCaD(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray_extCaD(:,1) .* plotArray_extCaD(:,1);
Area_values_normalized_noCa = Area_values_mutant ./ Area_values_mutant(1);

plot(t(1:indexBurstMutants(4)), Area_values_normalized_noCa(1:indexBurstMutants(4)), 'Color', mycolors(4,:), 'LineWidth', 2);
plot(t(indexBurstMutants(4):end), Area_values_normalized_noCa(indexBurstMutants(4):end), '--', 'Color', mycolors(4,:), 'LineWidth', 2);
plot(t(indexBurstMutants(4)), Area_values_normalized_noCa(indexBurstMutants(4)), 'o', 'Color', mycolors(4,:), 'LineWidth', 2);
plot([t(1) t(end)], [Area_values_normalized_noCa(indexBurstMutants(4)) Area_values_normalized_noCa(indexBurstMutants(4))], 'Color', mycolors(4,:), 'LineWidth', 1);

% tricalbin overexpression
Area_values_mutant = 2.0*pi^2*small_r.*plotArray_tricalbO(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray_tricalbO(:,1) .* plotArray_tricalbO(:,1);
Area_values_normalized_tricalbO = Area_values_mutant ./ Area_values_mutant(1);

plot(t(1:indexBurstMutants(5)), Area_values_normalized_tricalbO(1:indexBurstMutants(5)), 'Color', mycolors(5,:), 'LineWidth', 2);
plot(t(indexBurstMutants(5):end), Area_values_normalized_tricalbO(indexBurstMutants(5):end), '--', 'Color', mycolors(5,:), 'LineWidth', 2);
%plot(t(indexBurstMutants(5)), Area_values_normalized_mutant(indexBurstMutants(5)), 'o', 'Color', mycolors(5,:));


% no endocyotsis
Area_values_mutant = 2.0*pi^2*small_r.*plotArray_endoD(:,1)  + 4*pi*small_r^2 + 2*pi.*plotArray_endoD(:,1) .* plotArray_endoD(:,1);
Area_values_normalized_endoD = Area_values_mutant ./ Area_values_mutant(1);

plot(t(1:indexBurstMutants(6)), Area_values_normalized_endoD(1:indexBurstMutants(6)), 'Color', mycolors(6,:), 'LineWidth', 2);
plot(t(indexBurstMutants(6):end), Area_values_normalized_endoD(indexBurstMutants(6):end), '--', 'Color', mycolors(6,:), 'LineWidth', 2);
plot(t(indexBurstMutants(6)), Area_values_normalized_endoD(indexBurstMutants(6)), 'o', 'Color', mycolors(6,:), 'LineWidth', 2);
plot([t(1) t(end)], [Area_values_normalized_endoD(indexBurstMutants(6)) Area_values_normalized_endoD(indexBurstMutants(6))], 'Color', mycolors(6,:), 'LineWidth', 1);


%set(findall(gca, 'Type', 'Line'),'LineWidth',2);

ylabel('Relative surface area', 'FontSize', 16)

xlabel('Time (s)', 'FontSize', 20)
box on

figure(5)
clf()
hold on

x_values = t;
y_values = [1:length(indexBurstMutants)+1];


[X, Y] = meshgrid(x_values, y_values);

%color_values = [plotArray(:,12), plotArray_pil1D(:,12), plotArray_tricalbD(:,12), plotArray_exoD(:,12), plotArray_extCaD(:,12), plotArray_tricalbO(:,12), plotArray_endoD(:,12)];
color_values = [plotArray_endoD(:,12), plotArray_tricalbO(:,12), plotArray_extCaD(:,12), plotArray_exoD(:,12), plotArray_tricalbD(:,12), plotArray_pil1D(:,12), plotArray(:,12)];
color_values = color_values';

ax = gca();
ax.FontSize = 14;

%{
disp(size(X))
disp(size(Y))
disp(size(color_values))

disp(plotArray_pil1D(:,12));
%}

%h = pcolor(X, Y, color_values);
%set(h, 'EdgeColor', 'none');

%y_values = ["WT", "pil1D", "tricalbD", "exoD", "extCaD", "tricalb0", "endoD"];

%{
img = imagesc(x_values,y_values,color_values);
colormap(flipud('gray'));
colorbar()
val = 0.000001;

B = bwboundaries(img);

clim([epsThres-val, epsThres+val])

%}

%colormap('parula')
colormap('jet')

imagesc(x_values,y_values,color_values);

%disp(indexThreshold_tricalbO)

%draw rectangles
for i = 1:2:length(indexThreshold)
    if length(indexThreshold) - i == 0
        if indexThreshold(1) ~= -1
            rectangle('Position', [indexThreshold(i)*0.01 6.5 (length(t)-indexThreshold(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold(i)*0.01 6.5 (indexThreshold(i+1)-indexThreshold(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end
for i = 1:2:length(indexThreshold_pil1D)
    if length(indexThreshold_pil1D) - i == 0
        if indexThreshold_pil1D(1) ~= -1
            rectangle('Position', [indexThreshold_pil1D(i)*0.01 5.5 (length(t)-indexThreshold_pil1D(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold_pil1D(i)*0.01 5.5 (indexThreshold_pil1D(i+1)-indexThreshold_pil1D(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end
for i = 1:2:length(indexThreshold_tricalbD)
    if length(indexThreshold_tricalbD) - i == 0
        if indexThreshold_tricalbD(1) ~= -1
            rectangle('Position', [indexThreshold_tricalbD(i)*0.01 4.5 (length(t)-indexThreshold_tricalbD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold_tricalbD(i)*0.01 4.5 (indexThreshold_tricalbD(i+1)-indexThreshold_tricalbD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end
for i = 1:2:length(indexThreshold_exoD)
    if length(indexThreshold_exoD) - i == 0
        if indexThreshold_exoD(1) ~= -1
            rectangle('Position', [indexThreshold_exoD(i)*0.01 3.5 (length(t)-indexThreshold_exoD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold_exoD(i)*0.01 3.5 (indexThreshold_exoD(i+1)-indexThreshold_exoD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end
for i = 1:2:length(indexThreshold_extCaD)
    if length(indexThreshold_extCaD) - i == 0
        if indexThreshold_extCaD(1) ~= -1
            rectangle('Position', [indexThreshold_extCaD(i)*0.01 2.5 (length(t)-indexThreshold_extCaD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold_extCaD(i)*0.01 2.5 (indexThreshold_extCaD(i+1)-indexThreshold_extCaD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end
for i = 1:2:length(indexThreshold_tricalbO)
    if length(indexThreshold_tricalbO) - i == 0
        if indexThreshold_tricalbO(1) ~= -1
            rectangle('Position', [indexThreshold_tricalbO(i)*0.01 1.5 (length(t)-indexThreshold_tricalbO(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold_tricalbO(i)*0.01 1.5 (indexThreshold_tricalbO(i+1)-indexThreshold_tricalbO(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end
for i = 1:2:length(indexThreshold_endoD)
    if length(indexThreshold_endoD) - i == 0
        if indexThreshold_endoD(1) ~= -1
            rectangle('Position', [indexThreshold_endoD(i)*0.01 0.5 (length(t)-indexThreshold_endoD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
        end
    else
        rectangle('Position', [indexThreshold_endoD(i)*0.01 0.5 (indexThreshold_endoD(i+1)-indexThreshold_endoD(i))*0.01 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
    end
end

%{
rectangle('Position', [indexThreshold(1)*0.01 0.5 (indexThreshold(2)-indexThreshold(1))*0.01 1.0], 'LineStyle', '--')
rectangle('Position', [indexThreshold_pil1D(1)*0.01 1.5 (indexThreshold_pil1D(2)-indexThreshold_pil1D(1))*0.01 1.0], 'LineStyle', '--')
rectangle('Position', [indexThreshold_tricalbD(1)*0.01 2.5 (indexThreshold_tricalbD(2)-indexThreshold_tricalbD(1))*0.01 1.0], 'LineStyle', '--')
rectangle('Position', [indexThreshold_exoD(1)*0.01 3.5 (indexThreshold_exoD(2)-indexThreshold_exoD(1))*0.01 1.0], 'LineStyle', '--')
rectangle('Position', [indexThreshold_extCaD(1)*0.01 4.5 (indexThreshold_extCaD(2)-indexThreshold_extCaD(1))*0.01 1.0], 'LineStyle', '--')
rectangle('Position', [indexThreshold_tricalbO(1)*0.01 5.5 (indexThreshold_tricalbO(2)-indexThreshold_tricalbO(1))*0.01 1.0], 'LineStyle', '--')
rectangle('Position', [indexThreshold_endoD(1)*0.01 6.5 (indexThreshold_endoD(2)-indexThreshold_endoD(1))*0.01 1.0], 'LineStyle', '--')
%}

cb = colorbar();
clim([0.0, 0.10]);
ylabel(cb, 'Strain', 'FontSize', 20)

xlim([0.0, 600]);
ylim([0.5, 7.5]);

xlabel("Time (s)", 'FontSize', 20)
ylabel("")

yticklabels(flip(["WT", "{\it pil1}\Delta", "{\ittcb1}\Delta{\ittcb2}\Delta", "BFA", "no external Ca", "{\ittcb1}{\ittcb2} OE", "LatA"]));
ax = gca;
ax.YAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

box on

figure(6)
clf()
hold on




% need to bin relative area values
x_values = 1.0:0.005:1.4;
y_values = [1:length(indexBurstMutants)+1];

%color_values = [plotArray_endoD(:,12), plotArray_tricalbO(:,12), plotArray_extCaD(:,12), plotArray_exoD(:,12), plotArray_tricalbD(:,12), plotArray_pil1D(:,12), plotArray(:,12)];
%color_values = [zeros(1,length(x_values)), zeros(1,length(x_values)), zeros(1,length(x_values)), zeros(1,length(x_values)), zeros(1,length(x_values)), zeros(1,length(x_values)), zeros(1,length(x_values))];

color_values = zeros(7, length(x_values));


for i = 1:length(x_values)-1
    count = 0;
    for j = 1:length(Area_values_normalized)
        if Area_values_normalized(j) >= x_values(i) && Area_values_normalized(j) <= x_values(i+1)
            color_values(7,i) = plotArray(j,12) + color_values(7,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(7,i) = color_values(7,i) / count;
    end

    %pil1D
    count = 0;
    for j = 1:length(Area_values_normalized_pil1D)
        if Area_values_normalized_pil1D(j) >= x_values(i) && Area_values_normalized_pil1D(j) <= x_values(i+1)
            color_values(6,i) = plotArray_pil1D(j,12) + color_values(6,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(6,i) = color_values(6,i) / count;
    end

    %tricalbD
    count = 0;
    for j = 1:length(Area_values_normalized_tricalbD)
        if Area_values_normalized_tricalbD(j) >= x_values(i) && Area_values_normalized_tricalbD(j) <= x_values(i+1)
            color_values(5,i) = plotArray_tricalbD(j,12) + color_values(5,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(5,i) = color_values(5,i) / count;
    end

    %exoD
    count = 0;
    for j = 1:length(Area_values_normalized_exoD)
        if Area_values_normalized_exoD(j) >= x_values(i) && Area_values_normalized_exoD(j) <= x_values(i+1)
            color_values(4,i) = plotArray_exoD(j,12) + color_values(4,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(4,i) = color_values(4,i) / count;
    end

    %extCaD
    count = 0;
    for j = 1:length(Area_values_normalized_noCa)
        if Area_values_normalized_noCa(j) >= x_values(i) && Area_values_normalized_noCa(j) <= x_values(i+1)
            color_values(3,i) = plotArray_extCaD(j,12) + color_values(3,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(3,i) = color_values(3,i) / count;
    end

    %tricalbOE
    count = 0;
    for j = 1:length(Area_values_normalized_tricalbO)
        if Area_values_normalized_tricalbO(j) >= x_values(i) && Area_values_normalized_tricalbO(j) <= x_values(i+1)
            color_values(2,i) = plotArray_tricalbO(j,12) + color_values(2,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(2,i) = color_values(2,i) / count;
    end

    %endoD
    count = 0;
    for j = 1:length(Area_values_normalized_endoD)
        if Area_values_normalized_endoD(j) >= x_values(i) && Area_values_normalized_endoD(j) <= x_values(i+1)
            color_values(1,i) = plotArray_endoD(j,12) + color_values(1,i);
            count = count + 1;
        end
    end
    if count > 0
        color_values(1,i) = color_values(1,i) / count;
    end

end









disp(size(color_values))
disp(size(x_values))
disp(size(y_values))
%color_values = color_values';


colormap('jet')

imagesc(x_values,y_values,color_values);

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

cb = colorbar();
clim([0.0, 0.10]);
ylabel(cb, 'Strain', 'FontSize', 20)

xlim([1.0, 1.4]);
ylim([0.5, 7.5]);

xlabel("Relative surface area", 'FontSize', 20)
ylabel("")

yticklabels(flip(["WT", "{\it pil1}\Delta", "{\ittcb1}\Delta{\ittcb2}\Delta", "BFA", "no external Ca", "{\ittcb1}{\ittcb2} OE", "LatA"]));
box on

%draw rectangles
for i = 1:size(color_values,1)
    starting_pos = -1;

    for j = 1:size(color_values, 2)
        if color_values(i, j) > epsThres && starting_pos == -1
            starting_pos = x_values(j);
        elseif color_values(i, j) < epsThres && starting_pos > -1
            % now draw the rectangle
            disp(i)
            rectangle('Position', [starting_pos i-0.5 x_values(j)-starting_pos 1.0], 'LineStyle', '--', 'LineWidth', 1.5)
            starting_pos = -1;
        end
    end
end