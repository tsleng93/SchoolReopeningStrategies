load('Exampleoutput');

tb = datetime(2020, 08, 24);
ts = datetime(2021, 5, 23);
t = tb:ts;



%{
load('Jul4currentstrategyworkspace');
load('Jul4noisolworkspace');
load('Jul4nomasstestingworkspace');
load('Jul4SCTworkspace');
%}

AAA = [];





col1 = [0.00, 0.45, 0.74];
col2 = '#ff5765';
col3 = '#228b22';
col4 = '#AF58BA';


ciup = 0.975;
cidown = 0.025;






%Figure 3
figure;
A = Rinfs2_vec1./Rinfs1_vec1;
above = quantile((A), ciup);
below = quantile(A, cidown);

above = movmean(above, 7);
below = movmean(below, 7);


plot(t(8:end), mean(A), 'color', '[0.79 0.35 0.83]'); hold on
p = patch([0:265 fliplr(0:265)], [below above(end:-1:1)], 'r', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.79 0.35 0.83];


AA = movmean(mean(A), 7);
plot(t(8:end), AA, 'color', '#AF58BA',  'LineWidth', 1.5); hold on
%p = patch([0:265 fliplr(0:265)], [below above(end:-1:1)], r, 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
%p.FaceColor = '#AF58BA';

xlabel('Date');
ylabel('R_{school}');

%legend({'no mass testing, isolating close contacts', 'mass testing, isolating close contacts', 'mass testing, isolating year groups', 'R_{school} = 1', 'Schools closed'});
plot([t(1) t(end)], [1 1], 'k--');
 set(gcf, 'Position', [300, 300, 600, 400]);
 set(gca, 'fontsize', 12);
 
 h=fill([t(17*7+1),t(27*7+7),t(27*7+7),t(17*7+1)],[0,0,1.4,1.4],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];


h=fill([t(31*7+1),t(33*7+7),t(33*7+7),t(31*7+1)],[0,0, 1.4,1.4],'red');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

h=fill([t(9*7+1),t(9*7+7),t(9*7+7),t(9*7+1)],[0,0,1.4,1.4],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

xlim([t(8) t(end)]);
legend({'Model (daily)', 'Model (weekly moving average)', 'R_{school} = 1', 'School closures/holidays'});



figure;

A = InternalInc_vec1 + ExternalInc_vec1;
% % of cases in schools over half-terms

AA = sum(InternalInc_vec1(:, 1:54), 2)./sum(A(:, 1:54), 2);

mean(AA)
quantile(AA, ciup)
quantile(AA, cidown)

AA = sum(InternalInc_vec1(:, 65:111), 2)./sum(A(:, 65:111), 2);
mean(AA)
quantile(AA, ciup)
quantile(AA, cidown)

AA = sum(InternalInc_vec1(:, 191:265), 2)./sum(A(:, 191:265), 2);
mean(AA)
quantile(AA, ciup)
quantile(AA, cidown)

 A = InternalInc_vec1 + ExternalInc_vec1;
above = 100*quantile((A), ciup);
below = 100*quantile(A, cidown);
 
AA = 100*mean(A);
plot(t(8:end), AA, 'color', '#009ADE',  'LineWidth', 1); hold on
p = patch([0:265 fliplr(0:265)], [below above(end:-1:1)], 'r', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = '#009ADE';

 A = ExternalInc_vec1;
above = 100*quantile((A), ciup);
below = 100*quantile(A, cidown);
 
AA = 100*mean(A);
plot(t(8:end), AA, 'color', '#FF1F5B',  'LineWidth', 1); hold on
p = patch([0:265 fliplr(0:265)], [below above(end:-1:1)], 'r', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = '#FF1F5B';


h=fill([t(17*7+1),t(27*7+7),t(27*7+7),t(17*7+1)],[0,0,1.4,1.4],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];


h=fill([t(31*7+1),t(33*7+7),t(33*7+7),t(31*7+1)],[0,0, 1.4,1.4],'red');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

h=fill([t(9*7+1),t(9*7+7),t(9*7+7),t(9*7+1)],[0,0,1.4,1.4],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

xlim([t(8) t(end)]);

xlabel('Date');
ylabel('Incidence (%)');

ylim([0 0.3]);

legend({'All infections', 'External infections','Schools closures/holidays'});
 set(gcf, 'Position', [300, 300, 600, 400]);
 set(gca, 'fontsize', 12);



 
 
 %Figure 2
 figure;
 plot(t(8:end), 100*mean(tests_vec1), 'color', col4, 'LineWidth', 2);
A = tests_vec1;
above = quantile((A), ciup);
below = quantile(A, cidown);
p = patch([0:265 fliplr(0:265)], 100*[below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = '#AF58BA';

h=fill([t(17*7+1),t(27*7+7),t(27*7+7),t(17*7+1)],[0,0,0.09,0.09],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];


h=fill([t(31*7+1),t(33*7+7),t(33*7+7),t(31*7+1)],[0,0, 0.09,0.09],'red');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

h=fill([t(9*7+1),t(9*7+7),t(9*7+7),t(9*7+1)],[0,0,0.09,0.09],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

xlim([t(8) t(end)]);
 legend({'Model','Schools closures/holidays'});
 set(gcf, 'Position', [300, 300, 600, 400]);
 xlabel('Date');
 set(gca, 'fontsize', 12);
 
 
 figure;
  plot(t(end-76:end), 100*mean(poslft_vec1(:, (end-76):end)), 'color', col4, 'LineWidth', 2);

A = poslft_vec1(:, (end-76):end);

above = quantile((A), ciup);
below = quantile(A, cidown);
p = patch([0:76 fliplr(0:76)], 100*[below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = '#AF58BA';

h=fill([t(17*7+1),t(27*7+7),t(27*7+7),t(17*7+1)],[0,0,0.035,0.035],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];


h=fill([t(31*7+1),t(33*7+7),t(33*7+7),t(31*7+1)],[0,0, 0.035,0.035],'red');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

h=fill([t(9*7+1),t(9*7+7),t(9*7+7),t(9*7+1)],[0,0,0.035,0.035],'red', 'HandleVisibility', 'off');
h.FaceAlpha = 0.25;
h.EdgeColor = 'none';
h.FaceColor = [0.65 0.65 0.65];

xlim([t(end-76) t(end)]);
 legend({'Model','Schools closures/holidays'});

 set(gcf, 'Position', [300, 300, 600, 400]);
 xlabel('Date');
 ylabel('% pupils testing LFT positive');
 set(gca, 'fontsize', 12);
 
%}
 
 
 figure;
 x = 0:30;

above =  quantile(100*(hist_vec1./2979), ciup);
below =  quantile(100*(hist_vec1./2979), cidown);
above(31) = sum(above(31:end));
above = above(1:31);
below(31) = sum(below(31:end));
below = below(1:31);
y = [below; above];
hB=bar(0:30, [y(1,:); diff(y)].', ...
 'stacked','FaceColor','flat','EdgeColor','w'); hold on
hB(1).CData=ones(length(x),3);  % set the color of bottom to background (white)
%hB(2).CData = '#AF58BA';
hB(1).HandleVisibility = 'off';
hB(2).FaceColor = '#AF58BA';
%plot(0:30, 100*pillar_data_hist_1/sum(pillar_data_hist_1), 'k.', 'MarkerSize', 20);

 xticks([0 5 10 15 20 25 30]);
 xticklabels({'0', '5', '10', '15', '20', '25', '30+'});

 set(gcf, 'Position', [300, 300, 600, 400]);
 xlabel('Peak confirmed cases, N');
 ylabel('% of schools with peak N cases');
  legend({'Model (Sep-Dec)', 'Data (Sep-Dec)'});

 set(gca, 'fontsize', 12);
 %}
 
 
  figure;
 x = 0:30;

above =  quantile(100*(hist_vec2./2979), ciup);
below =  quantile(100*(hist_vec2./2979), cidown);
above(11) = sum(above(11:end));
above = above(1:11);
below(11) = sum(below(11:end));
below = below(1:11);
y = [below; above];
hB=bar(0:10, [y(1,:); diff(y)].', ...
 'stacked','FaceColor','flat','EdgeColor','w'); hold on
hB(1).CData=ones(length(x),3);  % set the color of bottom to background (white)
%hB(2).CData = '#AF58BA';
hB(1).HandleVisibility = 'off';
hB(2).FaceColor = '#AF58BA';
%plot(0:10, 100*pillar_data_hist_2/sum(pillar_data_hist_2), 'k.', 'MarkerSize', 20);

 set(gcf, 'Position', [300, 300, 600, 400]);
 
 xticks([0 2 4 6 8 10]);
 xticklabels({'0', '2', '4', '6', '8', '10+'});
 xlabel('Peak confirmed cases, N');
 ylabel('% of schools with peak N cases');
  legend({'Model (Mar-May)', 'Data (Mar-May)'});

 set(gca, 'fontsize', 12);
 %}

