% Code file for Figures 3.3 and 3.4
% Visualizes random and polar grids

clear, clc, close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultlegendinterpreter','latex');  

% Enable the use of customized functions generating the grid points
addpath('grids/');

%% Setup

% Switch flag for saving results to txt-file
save_results = 0;

% Set grid size
R = 12;
T = 2*R;

%% Generate grids

% Random grid
tic
x_rand = rand(R^2,2)-1/2;
t_rand = toc;
            
% Jittered grid
tic
x_jit = (-0.5+1/(2*R):1/R:0.5-1/(2*R));
[X1,X2] = ndgrid(x_jit,x_jit);
x_jit = [X1(:) X2(:)];
x_jit = x_jit + 1/(4*R)*(2*rand(R^2,2)-1);
clear X1 X2;
t_jit = toc;
            
% Polar grid
tic
x_p = polar(R,T);
t_polar = toc;
x_p = unique(x_p,'stable','rows');

% Modified polar grid
tic
x_mp = mpolar(R,T);
t_mpolar = toc;
x_mp = unique(x_mp,'stable','rows');

% Linogram grid
tic
x_l = linogram(R,T);
t_linogram = toc;
x_l = unique(x_l,'stable','rows');

% Golden angle polar grid
tic
x_gap = golden_angle_polar(R,T);
t_gap = toc;
x_gap = unique(x_gap,'stable','rows');

% Golden angle linogram grid
tic
x_gal = golden_angle_linogram(R,T);
t_gal = toc;
x_gal = unique(x_gal,'stable','rows');

% Spiral grid
tic
x_spir = spiral(R*T);
t_spir = toc;
x_spir = unique(x_spir,'stable','rows');

%% Visualization

figure(1); subplot(1,2,1); plot(x_rand(:,1),x_rand(:,2),'o'); title('Random grid (3.86)');
subplot(1,2,2); plot(x_jit(:,1),x_jit(:,2),'o'); title('Jittered grid (3.87)');
sgtitle('Figure 3.3: Exemplary grids of random kind of size~$N_1=N_2=12$.');

figure(2); subplot(2,3,1); plot(x_p(:,1),x_p(:,2),'o'); title('Polar grid (3.88)');
subplot(2,3,2); plot(x_mp(:,1),x_mp(:,2),'o'); title('Modified polar grid (3.89)');
subplot(2,3,3); plot(x_l(:,1),x_l(:,2),'o'); title('Linogram/pseudo-polar grid (3.90)');
subplot(2,3,4); plot(x_gap(:,1),x_gap(:,2),'o'); title('Golden angle polar grid, cf. (3.91)');
subplot(2,3,5); plot(x_gal(:,1),x_gal(:,2),'o'); title('Golden angle linogram grid (3.92)');
subplot(2,3,6); plot(x_spir(:,1),x_spir(:,2),'o'); title('Spiral grid (3.93)');
sgtitle('Figure 3.4: Exemplary grids of polar kind of size~\mbox{$R=12$} and~\mbox{$T=2R$}.');

%% Generate tables for tikz

if (save_results == 1)
fileID = fopen('grids.txt','w');
format = '%g %g \n';

fprintf(fileID,'Polar\n\n'); 
matrix = [x_p(:,1),x_p(:,2)];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n');

fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,'Modified polar\n\n'); 
matrix = [x_mp(:,1),x_mp(:,2)];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,'Linogram\n\n'); 
matrix = [x_l(:,1),x_l(:,2)];
fprintf(fileID,format,transpose(matrix));
fprintf(fileID,'\n\n');

fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,'Golden angle polar\n\n'); 
matrix = [x_gap(:,1),x_gap(:,2)];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,'Golden angle linogram\n\n'); 
matrix = [x_gal(:,1),x_gal(:,2)];
fprintf(fileID,format,transpose(matrix));

fprintf(fileID,'--------------------------------------------------------\n'); 
fprintf(fileID,'Spiral\n\n'); 
matrix = [x_spir(:,1),x_spir(:,2)];
fprintf(fileID,format,transpose(matrix));

fclose(fileID);
end%if