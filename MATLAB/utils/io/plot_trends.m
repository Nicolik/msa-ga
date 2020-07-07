function plot_trends(stats, id, description)
best_value = stats.best_value;
worst_value = min(stats.worst_values);
std_avg = std(stats.avg_values);
fprintf("Best value = %.1f\n", best_value);

figure('units','normalized','outerposition',[0 0 1 1])
xmin = 1;
xmax = size(stats.best_values,1);
plot(xmin:xmax, stats.best_values,  '-o', ...
     xmin:xmax, stats.avg_values,   ':*', ...
     xmin:xmax, stats.worst_values, ':d', ...
     'LineWidth', 1.5)
ymin = worst_value - 2 * std_avg;
ymax = best_value + 2 * std_avg;
ylim([ymin, ymax])
xtext = xmin+(xmax-xmin)*0.7;
ytext = ymin+(ymax-ymin)*0.2;
title("Optimization Process Trends")
xlabel("Generation")
ylabel("Fitness")
legend("Best Values", "Average Values", "Worst Values")
text(xtext,ytext,description,'FontName','FixedWidth')
txt_best = strcat({'\downarrow '}, sprintf('F: %.2f @ G: %03d',...
                                    best_value, stats.best_gen));
xtextb = stats.best_gen - (xmax-xmin) * 0.001;
ytextb = stats.best_value + (ymax-ymin) * 0.02; 
text(xtextb,ytextb,txt_best)
saveas(gcf,sprintf('./MATLAB/images/trends_%d_C%d.png',id,stats.chromosomes))
end

