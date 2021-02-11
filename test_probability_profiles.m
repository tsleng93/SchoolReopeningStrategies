%This .m file produces the figures for the test probability profiles of PCR
%and LFTs (Supplementary Figure S1)


%symptomatic testing profile
%Lateral flow
LatTable = readtable('lat_Curve_summary.csv');
lat_test_sym = table2array(LatTable(:,2:4));

%PCR
PCRTable = readtable('PCR_Curve_summary.csv');
PCR_test_sym = table2array(PCRTable(:,2:4));

%times
Symx = table2array(PCRTable(:,1));

%asymptomatic
Asymx = Symx;
Peak = find(PCR_test_sym == max(PCR_test_sym),1);
Asymx((Peak+1):end) = Asymx(Peak) + (6.7/10.5)*(Asymx((Peak+1):end) - Asymx(Peak));

%Find values by interpolating for asymptomatic
N = (10*Asymx(end) +1);

lat_test_asym(1:N,1) = interp1(Asymx, lat_test_sym(:,1), 0.1*(0:(N-1)));
lat_test_asym(1:N,2) = interp1(Asymx, lat_test_sym(:,2), 0.1*(0:(N-1)));
lat_test_asym(1:N,3) = interp1(Asymx, lat_test_sym(:,3), 0.1*(0:(N-1)));

%fit final values
fitlf1 = fit(Asymx(201:end),lat_test_sym(201:end, 1), 'Gauss1');
lat_test_asym((N+1):301, 1) = fitlf1(0.1*(N:300));
fitlf2 = fit(Asymx(201:end),lat_test_sym(201:end, 2), 'Gauss1');
lat_test_asym((N+1):301, 2) = fitlf2(0.1*(N:300));
fitlf3 = fit(Asymx(201:end),lat_test_sym(201:end, 3), 'Gauss1');
lat_test_asym((N+1):301, 3) = fitlf3(0.1*(N:300));

PCR_test_asym(1:N,1) = interp1(Asymx, PCR_test_sym(:,1), 0.1*(0:(N-1)));
PCR_test_asym(1:N,2) = interp1(Asymx, PCR_test_sym(:,2), 0.1*(0:(N-1)));
PCR_test_asym(1:N,3) = interp1(Asymx, PCR_test_sym(:,3), 0.1*(0:(N-1)));

%fit final values
fitPCR1 = fit(Asymx(201:end),PCR_test_sym(201:end, 1), 'Gauss1');
PCR_test_asym(N:301, 1) = fitPCR1(0.1*(N-1:300));
fitPCR2 = fit(Asymx(201:end),PCR_test_sym(201:end, 2), 'Gauss1');
PCR_test_asym(N:301, 2) = fitPCR2(0.1*((N-1):300));
fitPCR3 = fit(Asymx(201:end),PCR_test_sym(201:end, 3), 'Gauss1');
PCR_test_asym(N:301, 3) = fitPCR3(0.1*((N-1):300));


%write files
csvwrite('lat_Curve_asym.csv', lat_test_asym);
csvwrite('PCR_Curve_asym.csv', PCR_test_asym);



%Adjust for specificity
lat_test_sym(lat_test_sym(:,1) < 0.003, 1) = 0.003;
lat_test_sym(lat_test_sym(:,2) < 0.003, 2) = 0.003;
lat_test_sym(lat_test_sym(:,3) < 0.003, 3) = 0.003;
lat_test_asym(lat_test_asym(:,1) < 0.003, 1) = 0.003;
lat_test_asym(lat_test_asym(:,2) < 0.003, 2) = 0.003;
lat_test_asym(lat_test_asym(:,3) < 0.003, 3) = 0.003;


figure; %Fig S1B
set(gcf,'Position',[300,300,600,400] );
plot(Symx, PCR_test_sym(:,2), 'color', [0.06, 1, 1], 'linewidth', 1); hold on
plot(Symx, PCR_test_asym(:,2), '--', 'color', [0.06, 1, 1], 'linewidth', 1); 
plot(Symx, PCR_test_sym(:,1), 'color', [0.30 0.75, 0.93], 'linewidth', 1); hold on
plot(Symx, PCR_test_asym(:,1), '--', 'color', [0.30 0.75 0.93], 'linewidth', 1);
plot(Symx, PCR_test_sym(:,3), 'color', [0.00 0.45 0.74], 'linewidth', 1); hold on
plot(Symx, PCR_test_asym(:,3), '--', 'color', [0.00 0.45 0.74], 'linewidth', 1); 


xlabel('Days since infection');
ylabel('Probability of +ve PCR test');
ylim([0 0.9]);
legend('Symptomatic - low', 'Asymptomatic - low', 'Symptomatic - baseline', 'Asymptomatic - baseline', 'Symptomatic - high', 'Asymptomatic - high');
set(gca, 'fontsize', 12);

figure; %Figure S1A
set(gcf,'Position',[300,300,600,400] );

plot(Symx, lat_test_sym(:,2), 'color', [0.06, 1, 1], 'linewidth', 1); hold on
plot(Symx, lat_test_asym(:,2), '--', 'color', [0.06, 1, 1], 'linewidth', 1); 
plot(Symx, lat_test_sym(:,1), 'color', [0.30 0.75, 0.93], 'linewidth', 1); hold on
plot(Symx, lat_test_asym(:,1), '--', 'color', [0.30 0.75 0.93], 'linewidth', 1);
plot(Symx, lat_test_sym(:,3), 'color', [0.00 0.45 0.74], 'linewidth', 1); hold on
plot(Symx, lat_test_asym(:,3), '--', 'color', [0.00 0.45 0.74], 'linewidth', 1); 

xlabel('Days since infection');
ylabel('Probability of +ve LF test');
ylim([0 0.9]);
legend('Symptomatic - low', 'Asymptomatic - low', 'Symptomatic - baseline', 'Asymptomatic - baseline', 'Symptomatic - high', 'Asymptomatic - high');
set(gca, 'fontsize', 12);
