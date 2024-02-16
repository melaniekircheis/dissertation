% Code file for Figure 4.6

clear, clc, close all
fprintf('Started %s\n', datestr(datetime('now')))

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

%% Read tabulated data from [Medhurst,Roberts 65]

v = [2	1
4	0.666666666666667
6	0.55
8	0.479365079365079
10	0.430417768959436
12	0.393925565175565
14	0.365370869485453
16	0.342240261355337
18	0.323009394156998
20	0.306693101737992
22	0.292622687231527
24	0.280326198550056
26	0.269459771241317
28	0.259766148032194
30	0.2510485150
32	0.2431533907
34	0.2359590855
36	0.2293677077
38	0.2232994960
40	0.2176887196
42	0.2124806563
44	0.2076293306
46	0.2030957944
48	0.1988468037
50	0.1948537880
52	0.1910920377
54	0.1875400585
56	0.1841790530
58	0.1809925000
60	0.1779658125
62	0.1750860560
64	0.1723417150
66	0.1697224991
68	0.1672191807
70	0.1648234572
72	0.1625278358
74	0.1603255342
76	0.1582103967
78	0.1561768214
80	0.1542196975
82	0.1523343521
84	0.1505165022
86	0.1487622145
88	0.1470678691
90	0.1454301285
92	0.1438459096
94	0.1423123593
96	0.1408268331
98	0.1393868757
100	0.1379902041];

% Compute normalized vector
b0 = sqrt(v(:,1)).*v(:,2);

%% Visualization

figure(3); plot(2:50,b0(2:end),'b-',2:50,4/3*ones(49,1),'k-.',2:50,sqrt(6/pi)*ones(49,1),'k-.'); 
lgd = legend('$\sqrt{2s}\,B_{2s}(0)$'); xlabel('$s$'); lgd.FontSize = 11; 
yticks([4/3,sqrt(6/pi)]); yticklabels({'$\frac{4}{3}$','$\sqrt{\frac{6}{\pi}}$'});
xlim([2,50]); xticks([2,10,20,30,40,50]);
title(['Figure 4.6: The sequence~\mbox{$\big(\sqrt{2s}\,B_{2s}(0)\big)_{s=2}^{50}$}.'])

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('assertion_bspline.txt','w');
format = '%d %1.4e \n';
fprintf(fileID,'B-Spline at x=0\n\n');
matrix = [(2:50)',b0(2:end)];
fprintf(fileID,format,matrix.');
fprintf(fileID,'----------------------------------------\n');
matrix = [[2;50],4/3*ones(2,1)];
fprintf(fileID,format,matrix.');
fprintf(fileID,'----------------------------------------\n');
matrix = [[2;50],sqrt(6/pi)*ones(2,1)];
fprintf(fileID,format,matrix.');
fclose(fileID);
end%if