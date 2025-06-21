% units of molarity is M
delMolarity = 0.01;
delTimeDelay = 1.0;
molarities = 0.8:delMolarity:1.0;

disp(molarities)

% time delay in seconds
%timeDelaysExo = [0.0 5.0 10.0  15.0 30.0 45.0 60.0 90.0 120.0];
timeDelaysExo = 0.0:delTimeDelay:20;

k5_step_maxes = 0.0:1.0:2.0;
color_array = ["red", "blue", "green"];

[X,Y] = meshgrid(timeDelaysExo, 1.2-molarities);
%R = sqrt(X.^2 + Y.^2) + eps;
Z = zeros(length(molarities), length(timeDelaysExo));

baseStartingVolume = 9.6206E-17;

initialMolarityDelta = 1.0;
intialTimeDelay = 5.0;
initialk5Max = [1.5*2.0];
initialfACyto = [0.3];
initialfEis = [0.05];
initialtensionSValue3 = [0.01];
initialCaOut = 2.5E-3*1000;
colorList = ["blue" "red"];

%gap_height = 4.5e-6*1E6;
%small_r = gap_height*0.5;

hold off

epsilon_threshold = 0.056;
threshold_reached = zeros(length(molarities),1);


figure(1)
clf()

box on
Z = zeros(length(molarities), length(timeDelaysExo));
for mIndex = 1:length(molarities)
    for tIndex = 1:length(timeDelaysExo)
        %maxEpsValue = AreaBalancePassValues_20240321(molarities(mIndex), timeDelaysExo(tIndex), k5_step_maxes(kIndex));
        [maxEpsValue, plotArray, t] = AreaBalancePassValues_Final_20250325(molarities(mIndex), timeDelaysExo(tIndex), initialk5Max(1), baseStartingVolume, initialfEis(1), initialtensionSValue3(1), initialCaOut, true, true, false);

         
        fprintf("%f %f %f\n", molarities(mIndex), timeDelaysExo(tIndex), maxEpsValue)
        Z(mIndex, tIndex) = maxEpsValue;

        if maxEpsValue > epsilon_threshold && threshold_reached(mIndex) == 0
            threshold_reached(mIndex) = tIndex;
        elseif maxEpsValue < epsilon_threshold && threshold_reached(mIndex) ~= 0
            disp('Max eps dropped below threshold after having been above it' + str(mIndex))
        end
    end
end

x_values = timeDelaysExo;
y_values = molarities;
color_values = Z;

colormap('jet')

imagesc(x_values,y_values,color_values);

set(gca,'YDir','normal') 
%mesh(X, Y, Z, 'EdgeColor', color_array(kIndex), 'FaceAlpha', 0.0);

hold on



%mesh(X, Y, Z);

x_points = [];
y_points = [];
disp(threshold_reached)

for index = 1:length(molarities)

    %{
    if index > 1 && threshold_reached(index-1) ~= threshold_reached(index) && threshold_reached(index) > 0 && threshold_reached(index-1) > 0
        %plot([x_values(threshold_reached(index-1))-0.5*delTimeDelay, x_values(threshold_reached(index))-0.5*delTimeDelay], [y_values(index-1)+delMolarity*0.5, y_values(index-1)+delMolarity*0.5], 'k--')
        x_points = [x_points [x_values(threshold_reached(index-1))-0.5*delTimeDelay, x_values(threshold_reached(index))-0.5*delTimeDelay]];
        y_points = [y_points [y_values(index-1)+delMolarity*0.5, y_values(index-1)+delMolarity*0.5]];
    end

    if threshold_reached(index) > 1
        %plot([x_values(threshold_reached(index))-0.5*delTimeDelay, x_values(threshold_reached(index))-0.5*delTimeDelay], [y_values(index)+delMolarity*0.5, y_values(index)-delMolarity*0.5], 'k--')
        x_points = [x_points [x_values(threshold_reached(index))-0.5*delTimeDelay, x_values(threshold_reached(index))-0.5*delTimeDelay]];
        y_points = [y_points [y_values(index)+delMolarity*0.5, y_values(index)-delMolarity*0.5]];
    end
    %}

    if index > 1 && threshold_reached(index) ~= threshold_reached(index-1)
        %then draw a horizontal line connecting
        x_min_index = threshold_reached(index-1);
        x_min = 0.0;
        x_max = 0.0;

        if x_min_index == 0
            x_min = x_values(end)+0.5*delTimeDelay;
        elseif x_min_index == 1
            x_min = x_values(x_min_index)-0.5*delTimeDelay;
        else
            x_min = x_values(x_min_index)-0.5*delTimeDelay;
        end

        x_max_index = threshold_reached(index);
        if x_max_index == 0
            x_max = x_values(end)+0.5*delTimeDelay;
        elseif x_max_index == 1
            x_max = x_values(x_max_index)-0.5*delTimeDelay;
        else
            x_max = x_values(x_max_index)-0.5*delTimeDelay;
        end

        x_points = [x_points [x_min, x_max]];
        y_points = [y_points [y_values(index)-delMolarity*0.5, y_values(index)-delMolarity*0.5]];
    end

    if threshold_reached(index) > 1
        %then draw a straight vertical line
        x_points = [x_points [x_values(threshold_reached(index))-0.5*delTimeDelay, x_values(threshold_reached(index))-0.5*delTimeDelay]];
        y_points = [y_points [y_values(index)+delMolarity*0.5, y_values(index)-delMolarity*0.5]];
    end

    
end

disp(x_points)
disp(y_points)

plot(x_points, y_points, 'k--');

ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;

cb = colorbar();
clim([0.04, 0.085]);
ylabel(cb, 'Max strain', 'FontSize', 16)


%xlim([1.0, 1.4]);
%ylim([0.5, 7.5]);

ylabel("\Deltac (M)", 'FontSize', 20)
xlabel("Exocytosis time delay (s)")

%yticklabels(flip(["WT", "{\it pil1}\Delta", "{\ittcb1}\Delta{\ittcb2}\Delta", "BFA", "no external Ca", "{\ittcb1}{\ittcb2} OE", "LatA"]));
