%% Gonzalez, Bavassi & Laje (2019) "Response to perturbations as a built-in feature in a mathematical model for paced finger tapping"
% Submitted to PRE - October 2019
% Run this to plot all figures in the manuscript


%% Load behavioral data


clear all;

% load data
data_filename = 'perturb_exp.dat';
perturb_sizes_nobase = [50 40 30 20 10];
nbr_pert_sizes = length(perturb_sizes_nobase);
datos_mean_err = load(data_filename);	% columns: time,-50,-40,-30,-20,-10,10,20,30,40,50

% time series
datos_asyn = [];
datos_time = datos_mean_err(:,1);
datos_asyn(1,:,:) = datos_mean_err(:,[2:2:10]);		% negative steps
datos_asyn(2,:,:) = fliplr(datos_mean_err(:,[12:2:end]));	% positive steps
nbr_steps = size(datos_asyn,2);

% embedding
datos_asyn_embed = repmat(datos_asyn,[1 1 1 2]);
datos_asyn_embed(:,:,:,2) = datos_asyn_embed(:,:,:,2) - circshift(datos_asyn_embed(:,:,:,2),1,2);
datos_asyn_embed(:,1,:,:) = [];
postbl_ini = 23;
postbl_embed = mean(datos_asyn_embed(:,postbl_ini:end,:,:),2);
datos_asyn_embed_postbl = datos_asyn_embed - repmat(postbl_embed,[1 nbr_steps-1 1 1]);



%% Figure 1B: isochronous time series



load('asyn_iso_subj4_dis_size1:5_rep1:5.mat');
asyn_iso_ave = mean(asyn_iso);
asyn_iso_std = std(asyn_iso);
step_n = [1:length(asyn_iso)];


lwidth = 1;
msize = 15;
fsize = 12;
fsize_legend = 9;
fsize_text = 14;
fig_size_cm = [15 8];


figure(1);
clf(1);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

set(gca,'fontsize',fsize);
plot([0 step_n(end)],[0 0],'k-');
hold on;
plot(step_n,asyn_iso,'k.-','markersize',10);
set(gca,'xtick',[0:50:step_n(end)]);
set(gca,'ytick',[-40:20:40]);
xlim([0 step_n(end)]);
ylim([-45 45]);
xlabel('Step n');
ylabel('Asynchrony e_n (ms)');
text(-38,42,'(b)','fontsize',fsize_text);




filename = 'figure1b_isochronous';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 2: old experimental data


return_pts_x = [73.4679   30.9656;
60.1640   21.3633;
42.7626   16.8335;
27.9098    7.7074;
7.0368   -0.4513;
-13.2585   -3.1126;
-30.1763  -14.8603;
-46.4820  -19.2989;
-57.3800  -24.2211;
-65.7629  -32.1814];

return_pts_y = [28.0007  -27.1516;
30.1471  -21.5245;
23.5740  -14.3629;
-10.8972    2.0638;
-12.2941   13.9621;
-16.1259   22.0209;
-5.6086   34.2938;
-1.0791   44.5212];


lwidth = 1;
msize = 15;
fsize = 12;
fsize_legend = 9;
fsize_text = 14;
fig_size_cm = [15 25];

embed_ini = 11;
embed_fin = 20;
x_lims = [-5 20];
y_lims_neg = [-30 60];
y_lims_pos = [-60 40];
x_lims_embed = [-80 80];
y_lims_embed = [-60 60];

% new_cmap = ametrine(5);
new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);

figure(2);
clf(2);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

pert_sign = 1;
ptp = [];
subplot(4,1,1);
plot([0 0],y_lims_neg,'k-');
hold on;
plot(x_lims,[0 0],'k-');
% ptp = plot(datos_time,squeeze(datos_asyn(pert_sign,:,:)),'.--','linewidth',lwidth,'markersize',msize);
for pert_size = 1:nbr_pert_sizes
	ptp(pert_size) = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'.--','color',new_cmap(pert_size,:),'linewidth',lwidth,'markersize',msize);
end
xlim(x_lims);
ylim(y_lims_neg);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):5:x_lims(2)]);
set(gca,'ytick',[y_lims_neg(1):20:y_lims_neg(2)]+10);
xlabel('Step n relative to perturbation');
ylabel('e_n (ms)');
ptl = legend(ptp,'-50','-40','-30','-20','-10','location','northeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;
text(-8.5,60,'(a)','fontsize',fsize_text);


pert_sign = 2;
subplot(4,1,2);
plot([0 0],y_lims_pos,'k-');
hold on;
plot(x_lims,[0 0],'k-');
for pert_size = 1:nbr_pert_sizes
	ptp(pert_size) = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'.-','color',new_cmap(pert_size,:),'linewidth',lwidth,'markersize',msize);
end
xlim(x_lims);
ylim(y_lims_pos);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):5:x_lims(2)]);
set(gca,'ytick',[y_lims_pos(1):20:y_lims_pos(2)]);
xlabel('Step n relative to perturbation');
ylabel('e_n (ms)');
ptl = legend(ptp,'+50','+40','+30','+20','+10','location','southeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;
text(-8.5,40,'(b)','fontsize',fsize_text);

% embedding
subplot(4,1,[3 4]);
plot([0 0],y_lims_embed,'k-');
hold on;
plot(x_lims_embed,[0 0],'k-');
pert_sign = 1;
for pert_size = 1:nbr_pert_sizes
	plot(squeeze(datos_asyn_embed_postbl(pert_sign,embed_ini:embed_fin,pert_size,1)),squeeze(datos_asyn_embed_postbl(pert_sign,embed_ini:embed_fin,pert_size,2)),'.--','color',new_cmap(pert_size,:),'linewidth',lwidth,'markersize',msize);
end
hold on;
pert_sign = 2;
for pert_size = 1:nbr_pert_sizes
	plot(squeeze(datos_asyn_embed_postbl(pert_sign,embed_ini:embed_fin,pert_size,1)),squeeze(datos_asyn_embed_postbl(pert_sign,embed_ini:embed_fin,pert_size,2)),'.-','color',new_cmap(pert_size,:),'linewidth',lwidth,'markersize',msize);
end
ptp1 = plot(return_pts_x(:,1),return_pts_x(:,2),'ks-','linewidth',lwidth,'markersize',msize/4);
ptp2 = plot(return_pts_y(:,1),return_pts_y(:,2),'kd-','linewidth',lwidth,'markersize',msize/4);
set(ptp1,'MarkerFaceColor','k');
set(ptp2,'MarkerFaceColor','k');
xlim(x_lims_embed);
ylim(y_lims_embed);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims_embed(1):20:x_lims_embed(2)]);
set(gca,'ytick',[y_lims_embed(1):20:y_lims_embed(2)]);
xlabel('e_n (ms)');
ylabel('e_n - e_{n-1} (ms)');
text(-103,60,'(c)','fontsize',fsize_text);
legend([ptp1 ptp2],'return points (horiz)','return points (vert)','location','southeast');
legend boxoff;


filename = 'figure2_expdata';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 3: fitting results

% model parameters
a = 0.981;
b = 0.266;
c = -0.823;
d = 0.0238;
alpha = -2.213e-05;
beta = -7.845e-05;
gamma = 5.338e-05;
delta = 3.348e-3;
T0 = 500;		% ms


postbl = squeeze(mean(datos_asyn(:,postbl_ini:end,:),2));
exp_postbaseline = postbl(:,1)';

perturb_bip = 12;
y_ini = [0 T0 T0];


% postbaseline set to experimental values
perturb_sizes = [50];
nbr_perturb_sizes = length(perturb_sizes);
% variables
y = nan(2,nbr_perturb_sizes,nbr_steps,3);
asyn = nan(2,nbr_perturb_sizes,nbr_steps);
% parameter
T_all = T0*ones(2,nbr_perturb_sizes,nbr_steps);
T_all(1,:,perturb_bip:end) = T_all(1,:,perturb_bip:end) - repmat(perturb_sizes,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
T_all(2,:,perturb_bip:end) = T_all(2,:,perturb_bip:end) + repmat(perturb_sizes,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
for pert_sign = 1:2
	for pert_size = 1:nbr_perturb_sizes
		y(pert_sign,pert_size,1,:) = y_ini;
		perturb_size = perturb_sizes(pert_size);
		baseline = 0;
		for i = 1:nbr_steps-1
			p = y(pert_sign,pert_size,i,1);
			x = y(pert_sign,pert_size,i,2);
			s = y(pert_sign,pert_size,i,3);
			T = T_all(pert_sign,pert_size,i);
			if i==1
				Tprev = T0;
			else
				Tprev = T_all(pert_sign,pert_size,i-1);
			end
			if i==perturb_bip
				baseline = exp_postbaseline(pert_sign);
			end
			f = a*(p-(T-s)-baseline) + b*(x-T) ...
				+ alpha*(p-(T-s)-baseline)^3 ...
				+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
				+ gamma*(x-T)^3 ...
				+ baseline;
			g = c*(p-(T-s)-baseline) + d*(x-T) ...
				+ delta*(p-(T-s)-baseline)^2 ...
				+ T;
			h = T;
			
			DeltaT = T - Tprev;
			asyn(pert_sign,pert_size,i) = p - DeltaT;
			
			y(pert_sign,pert_size,i+1,1) = f;
			y(pert_sign,pert_size,i+1,2) = g;
			y(pert_sign,pert_size,i+1,3) = h;
		end
		asyn(pert_sign,pert_size,end) = p - DeltaT;
	end
end
y(y>1000) = NaN;
asyn(asyn>1000) = NaN;
x = y(:,:,:,2);
x_shift = x - T_all;



% postbaseline set to zero
perturb_sizes_nobase = [70 60 50 40 30 20 10];
nbr_perturb_sizes = length(perturb_sizes_nobase);
% variables
y = nan(2,nbr_perturb_sizes,nbr_steps,3);
asyn_nobase = nan(2,nbr_perturb_sizes,nbr_steps);
% parameter
T_all = T0*ones(2,nbr_perturb_sizes,nbr_steps);
T_all(1,:,perturb_bip:end) = T_all(1,:,perturb_bip:end) - repmat(perturb_sizes_nobase,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
T_all(2,:,perturb_bip:end) = T_all(2,:,perturb_bip:end) + repmat(perturb_sizes_nobase,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
for pert_sign = 1:2
	for pert_size = 1:nbr_perturb_sizes
		y(pert_sign,pert_size,1,:) = y_ini;
		perturb_size = perturb_sizes_nobase(pert_size);
		baseline = 0;
		for i = 1:nbr_steps-1
			p = y(pert_sign,pert_size,i,1);
			x = y(pert_sign,pert_size,i,2);
			s = y(pert_sign,pert_size,i,3);
			T = T_all(pert_sign,pert_size,i);
			if i==1
				Tprev = T0;
			else
				Tprev = T_all(pert_sign,pert_size,i-1);
			end
			f = a*(p-(T-s)-baseline) + b*(x-T) ...
				+ alpha*(p-(T-s)-baseline)^3 ...
				+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
				+ gamma*(x-T)^3 ...
				+ baseline;
			g = c*(p-(T-s)-baseline) + d*(x-T) ...
				+ delta*(p-(T-s)-baseline)^2 ...
				+ T;
			h = T;
			
			DeltaT = T - Tprev;
			asyn_nobase(pert_sign,pert_size,i) = p - DeltaT;
			
			y(pert_sign,pert_size,i+1,1) = f;
			y(pert_sign,pert_size,i+1,2) = g;
			y(pert_sign,pert_size,i+1,3) = h;
		end
		asyn_nobase(pert_sign,pert_size,end) = p - DeltaT;
	end
end
y(y>1000) = NaN;
asyn_nobase(asyn_nobase>1000) = NaN;
x_nobase = y(:,:,:,2);
x_nobase_shift = x_nobase - T_all;


lwidth = 1;
msize = 15;
lwidth2 = 2;
msize2 = 6;
fsize = 12;
fsize_legend = 10;
fsize_text = 14;
fig_size_cm = [15 20];

x_lims = [-5 20];
y_lims_neg = [-30 60];
y_lims_pos = [-60 40];
x_lims_embed = [-90 90];
y_lims_embed = [-90 90];
x_lims_eigen1 = [-90 39];
x_lims_eigen2 = [-90 25];
embed_ini = 12;
embed_fin = 20;

% eigenvectors
A = [[a b]; [c d]];
[M_eigenvec,M_eigenval] = eig(A);
lambda1 = M_eigenval(1,1);
lambda2 = M_eigenval(2,2);
eigenvec1 = M_eigenvec(:,1);
eigenvec2 = M_eigenvec(:,2);
slope1 = eigenvec1(2)/eigenvec1(1);
slope2 = eigenvec2(2)/eigenvec2(1);
vec1 = slope1*x_lims_eigen1;
vec2 = slope2*x_lims_eigen2;


new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);

figure(3);
clf(3);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

pert_size = 1;
subplot(9,1,[1 3]);
plot([0 0],[y_lims_pos(1) y_lims_neg(2)],'k-');
hold on;
plot(x_lims,[0 0],'k-');
pert_sign = 1;
cmap_color = 4;
ptp1 = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'o-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize2);
ptp2 = plot(datos_time,asyn(pert_sign,:),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
pert_sign = 2;
cmap_color = 1;
ptp3 = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'o-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize2);
ptp4 = plot(datos_time,asyn(pert_sign,:),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
xlim(x_lims);
ylim([y_lims_pos(1) y_lims_neg(2)]);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):5:x_lims(2)]);
set(gca,'ytick',[y_lims_pos(1):20:y_lims_neg(2)]);
xlabel('Step n relative to perturbation');
ylabel('e_n (ms)');
ptl = legend([ptp1 ptp2 ptp3 ptp4],'exp -50','fit -50','exp +50','fit +50','location','northeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;
text(-8.5,55,'(a)','fontsize',fsize_text);



msize2 = 10;
perturb_sizes_bound_pos_idx = [1 2 3 4 5 6 7];
perturb_sizes_bound_neg_idx = [1 2 3 4 5 6 7];
perturb_sizes_unbound_pos_idx = [];
perturb_sizes_unbound_neg_idx = [];
subplot(9,1,[5 9]);
plot([0 0],y_lims_embed,'k-');
hold on;
plot(x_lims_embed,[0 0],'k-');
ptp1 = plot(x_lims_eigen1,vec1,'k--');
plot(x_lims_eigen2,vec2,'k--');
pert_sign = 1;
cmap_color = 4;
for pert_size_idx = perturb_sizes_bound_neg_idx
	ptp2 = plot(squeeze(asyn_nobase(pert_sign,pert_size_idx,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size_idx,embed_ini:embed_fin)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize2);
end
pert_sign = 2;
cmap_color = 1;
for pert_size_idx = perturb_sizes_bound_pos_idx
	ptp3 = plot(squeeze(asyn_nobase(pert_sign,pert_size_idx,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size_idx,embed_ini:embed_fin)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize2);
end
xlim(x_lims_embed);
ylim(y_lims_embed);
set(gca,'fontsize',fsize);
set(gca,'xtick',[-80:40:80]);
set(gca,'ytick',[-80:40:80]);
xlabel('e_n (ms)');
ylabel('x_n - T_{post} (ms)');
ptl = legend([ptp2 ptp3 ptp1],'negative perturbations','positive perturbations','eigenvectors','location','southeast');
set(ptl,'fontsize',fsize_legend);%,'position',[0.27 0.13 0.3 0.05]);
legend boxoff;
text(-135,95,'(b)','fontsize',fsize_text);


filename = 'figure3_simuls';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);




%% Figure 4: types of phase space


lwidth = 1;
msize = 11;
lwidth2 = 2;
msize2 = 4;
fsize = 11;
fsize_legend = 10;
fsize_text = fsize;
fig_size_cm = [30 12];

x_lims = [-5 20];
y_lims_neg = [-30 60];
y_lims_pos = [-60 40];
x_lims_embed = [-90 90];
y_lims_embed = [-90 90];
embed_ini = 12;
embed_fin = 20;

new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);

figure(4);
clf(4);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);



for panel = 1:3

	% model parameters
	switch panel
		case 1
			a = 0.981;	% (loop=164)
			b = 0.266;
			c = -0.823;
			d = 0.0238;
			alpha = -2.213e-05;
			beta = -7.845e-05;
			gamma = 5.338e-05;
			delta = 3.348e-3;
			T0 = 500;		% ms
			% eigenvectors
			A = [[a b]; [c d]];
			[M_eigenvec,M_eigenval] = eig(A);
			lambda1 = M_eigenval(1,1);
			lambda2 = M_eigenval(2,2);
			eigenvec1 = M_eigenvec(:,1);
			eigenvec2 = M_eigenvec(:,2);
			slope1 = eigenvec1(2)/eigenvec1(1);
			slope2 = eigenvec2(2)/eigenvec2(1);
			vec1 = slope1*x_lims_embed;
			vec2 = slope2*x_lims_embed;
			panel_label = '(a)';
		case 2
			a = 0.8012;		% (loop=161)
			b = 0.6994;
			c = -0.22;
			d = 0.01668;
			alpha = -8.037e-05;
			beta = -8.058e-05;
			gamma = 9.764e-05;
			delta = 0.003485;
			T0 = 500;		% ms
			% eigenvectors
			A = [[a b]; [c d]];
			[M_eigenvec,M_eigenval] = eig(A);
			lambda1 = M_eigenval(1,1);
			lambda2 = M_eigenval(2,2);
			eigenvec1 = M_eigenvec(:,1);
			eigenvec2 = M_eigenvec(:,2);
			slope1 = eigenvec1(2)/eigenvec1(1);
			slope2 = eigenvec2(2)/eigenvec2(1);
			vec1 = slope1*x_lims_embed;
			vec2 = slope2*x_lims_embed;
			panel_label = '(b)';
		case 3
			a = 0.915;	% loop=3
			b = -0.055;
			c = 0.861;
			d = 0.364;
			alpha = 8.06e-05;
			beta = -9.22e-05;
			gamma = -3.484e-05;
			delta = -0.002776;
			T0 = 500;		% ms
			% eigenvectors
			A = [[a b]; [c d]];
			[M_eigenvec,M_eigenval] = eig(A);
			lambda1 = M_eigenval(1,1);
			lambda2 = M_eigenval(2,2);
			eigenvec1 = M_eigenvec(:,1);
			eigenvec2 = M_eigenvec(:,2);
			slope1 = eigenvec1(2)/eigenvec1(1);
			slope2 = eigenvec2(2)/eigenvec2(1);
			vec1 = slope1*x_lims_embed;
			vec2 = slope2*x_lims_embed;
			panel_label = '(c)';
	end

	postbl = squeeze(mean(datos_asyn(:,postbl_ini:end,:),2));
	exp_postbaseline = postbl(:,1)';

	perturb_bip = 12;
	y_ini = [0 T0 T0];


	% postbaseline set to experimental values
	perturb_sizes = [50];
	nbr_perturb_sizes = length(perturb_sizes);
	% variables
	y = nan(2,nbr_perturb_sizes,nbr_steps,3);
	asyn = nan(2,nbr_perturb_sizes,nbr_steps);
	% parameter
	T_all = T0*ones(2,nbr_perturb_sizes,nbr_steps);
	T_all(1,:,perturb_bip:end) = T_all(1,:,perturb_bip:end) - repmat(perturb_sizes,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
	T_all(2,:,perturb_bip:end) = T_all(2,:,perturb_bip:end) + repmat(perturb_sizes,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
	for pert_sign = 1:2
		for pert_size = 1:nbr_perturb_sizes
			y(pert_sign,pert_size,1,:) = y_ini;
			perturb_size = perturb_sizes(pert_size);
			baseline = 0;
			for i = 1:nbr_steps-1
				p = y(pert_sign,pert_size,i,1);
				x = y(pert_sign,pert_size,i,2);
				s = y(pert_sign,pert_size,i,3);
				T = T_all(pert_sign,pert_size,i);
				if i==1
					Tprev = T0;
				else
					Tprev = T_all(pert_sign,pert_size,i-1);
				end
				if i==perturb_bip
					baseline = exp_postbaseline(pert_sign);
				end
				f = a*(p-(T-s)-baseline) + b*(x-T) ...
					+ alpha*(p-(T-s)-baseline)^3 ...
					+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
					+ gamma*(x-T)^3 ...
					+ baseline;
				g = c*(p-(T-s)-baseline) + d*(x-T) ...
					+ delta*(p-(T-s)-baseline)^2 ...
					+ T;
				h = T;

				DeltaT = T - Tprev;
				asyn(pert_sign,pert_size,i) = p - DeltaT;

				y(pert_sign,pert_size,i+1,1) = f;
				y(pert_sign,pert_size,i+1,2) = g;
				y(pert_sign,pert_size,i+1,3) = h;
			end
			asyn(pert_sign,pert_size,end) = p - DeltaT;
		end
	end
	y(y>1000) = NaN;
	asyn(asyn>1000) = NaN;
	x = y(:,:,:,2);
	x_shift = x - T_all;



	% postbaseline set to zero
	perturb_sizes_nobase = [70 60 50 40 30 20 10];
	nbr_perturb_sizes = length(perturb_sizes_nobase);
	% variables
	y = nan(2,nbr_perturb_sizes,nbr_steps,3);
	asyn_nobase = nan(2,nbr_perturb_sizes,nbr_steps);
	% parameter
	T_all = T0*ones(2,nbr_perturb_sizes,nbr_steps);
	T_all(1,:,perturb_bip:end) = T_all(1,:,perturb_bip:end) - repmat(perturb_sizes_nobase,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
	T_all(2,:,perturb_bip:end) = T_all(2,:,perturb_bip:end) + repmat(perturb_sizes_nobase,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
	for pert_sign = 1:2
		for pert_size = 1:nbr_perturb_sizes
			y(pert_sign,pert_size,1,:) = y_ini;
			perturb_size = perturb_sizes_nobase(pert_size);
			baseline = 0;
			for i = 1:nbr_steps-1
				p = y(pert_sign,pert_size,i,1);
				x = y(pert_sign,pert_size,i,2);
				s = y(pert_sign,pert_size,i,3);
				T = T_all(pert_sign,pert_size,i);
				if i==1
					Tprev = T0;
				else
					Tprev = T_all(pert_sign,pert_size,i-1);
				end
				f = a*(p-(T-s)-baseline) + b*(x-T) ...
					+ alpha*(p-(T-s)-baseline)^3 ...
					+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
					+ gamma*(x-T)^3 ...
					+ baseline;
				g = c*(p-(T-s)-baseline) + d*(x-T) ...
					+ delta*(p-(T-s)-baseline)^2 ...
					+ T;
				h = T;

				DeltaT = T - Tprev;
				asyn_nobase(pert_sign,pert_size,i) = p - DeltaT;

				y(pert_sign,pert_size,i+1,1) = f;
				y(pert_sign,pert_size,i+1,2) = g;
				y(pert_sign,pert_size,i+1,3) = h;
			end
			asyn_nobase(pert_sign,pert_size,end) = p - DeltaT;
		end
	end
	y(y>1000) = NaN;
	asyn_nobase(asyn_nobase>1000) = NaN;
	x_nobase = y(:,:,:,2);
	x_nobase_shift = x_nobase - T_all;





	pert_size = 1;
	subplot(9,3,[panel panel+6]);
	plot([0 0],[y_lims_pos(1) y_lims_neg(2)],'k-');
	hold on;
	plot(x_lims,[0 0],'k-');
	pert_sign = 1;
	cmap_color = 4;
	ptp1 = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'o-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize2);
	plot(datos_time,asyn(pert_sign,:),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	pert_sign = 2;
	cmap_color = 1;
	ptp2 = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	plot(datos_time,asyn(pert_sign,:),'o-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize2);
	xlim(x_lims);
	ylim([y_lims_pos(1) y_lims_neg(2)]);
	set(gca,'fontsize',fsize);
	set(gca,'xtick',[x_lims(1):5:x_lims(2)]);
	set(gca,'ytick',[-50:25:50]);
	xlabel('Step n relative to perturbation');
	ylabel('e_n (ms)');
	if panel==1
		ptl = legend([ptp1 ptp2],'-50','+50','location','northeast');
		set(ptl,'fontsize',fsize_legend);
		legend boxoff;
		text(15,-38,'o  exp','fontsize',fsize_legend);
		text(15,-48,'\bullet  fit','fontsize',fsize_legend);
	end
	text(-11,60,panel_label,'fontsize',fsize_text);



	perturb_sizes_bound_pos_idx = [1 2 3 4 5 6];
	perturb_sizes_bound_neg_idx = [1 2 3 4 5 6];
	subplot(9,3,[panel+12 panel+24]);
	plot([0 0],y_lims_embed,'k-');
	hold on;
	plot(x_lims_embed,[0 0],'k-');
	ptp5 = plot(x_lims_embed,vec1,'k--');
	plot(x_lims_embed,vec2,'k--');
	pert_sign = 1;
	cmap_color = 4;
	for pert_size_idx = perturb_sizes_bound_neg_idx
		ptp1 = plot(squeeze(asyn_nobase(pert_sign,pert_size_idx,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size_idx,embed_ini:embed_fin)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	end
	pert_sign = 2;
	cmap_color = 1;
	for pert_size_idx = perturb_sizes_bound_pos_idx
		ptp3 = plot(squeeze(asyn_nobase(pert_sign,pert_size_idx,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size_idx,embed_ini:embed_fin)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	end
	xlim(x_lims_embed);
	ylim(y_lims_embed);
	set(gca,'fontsize',fsize);
	set(gca,'xtick',[-80:40:80]);
	set(gca,'ytick',[-80:40:80]);
	xlabel('e_n (ms)');
	ylabel('x_n - T_{post} (ms)');
	if panel==2
		ptl = legend([ptp1 ptp3 ptp5],'negative perturbs','positive perturbs','eigenvectors','location','southeast');
		set(ptl,'fontsize',fsize_legend);
		legend boxoff;
	end
end

filename = 'figure4_types';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 5: benefits

lwidth = 1;
msize = 8;
lwidth2 = 2;
msize2 = 5;
fsize = 11;
fsize_legend = 10;
fsize_text = fsize;
fig_size_cm = [30 12];


x_lims_step = [-5 20];
x_lims_sin_rand = [-5 30];
y_lims_step = [-20 60];
y_lims_step_period = [440 510];
y_lims_sin = [-100 100];
y_lims_sin_period = [450 550];
y_lims_rand = [-70 70];
y_lims_rand_period = [450 550];

new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);
data_color = 1;
delta_color = 4;

figure(5);
clf(5);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% model parameters
a = 0.981;
b = 0.266;
c = -0.823;
d = 0.0238;
alpha = -2.213e-05;
beta = -7.845e-05;
gamma = 5.338e-05;
delta = 3.348e-3;
T0 = 500;		% ms

perturb_bip = 12;
baseline = 0;
nbr_steps = 50;
datos_time_c = [1:nbr_steps]-perturb_bip;

for perturb = 1:3	%1: step change, 2: sinusoidal; 3: random

	y_ini = [0 T0 T0];
	% parameter
	T_all = T0*ones(1,nbr_steps);
	perturb_size_step = -50;
	perturb_size_sin = 47;%30;
	perturb_period = 15;
	perturb_size_rand = 30;
	if perturb==1
		T_all(perturb_bip:end) = T_all(perturb_bip:end) + perturb_size_step;
	elseif perturb==2
		T_all(perturb_bip:end) = T_all(perturb_bip:end) + perturb_size_sin*sin((2*pi/perturb_period)*([perturb_bip:nbr_steps]-perturb_bip+1));
	else
		T_all(perturb_bip:end) = T_all(perturb_bip:end) + perturb_size_rand*(2*rand(1,nbr_steps-perturb_bip+1)-1);
	end

	% variables
	y_predict = nan(nbr_steps,3);
	asyn_predict = nan(1,nbr_steps);
	y_byhand = nan(nbr_steps,3);
	asyn_byhand = nan(1,nbr_steps);
	e_minus = nan(1,nbr_steps);
	e_plus = nan(1,nbr_steps);
	y_naive = nan(nbr_steps,3);
	asyn_naive = nan(1,nbr_steps);

	y_predict(1,:) = y_ini;
	y_byhand(1,:) = y_ini;
	y_naive(1,:) = y_ini;
	for i = 1:nbr_steps-1
		T = T_all(i);
		if i==1
			Tprev = T0;
		else
			Tprev = T_all(i-1);
		end

		% our model
		p = y_predict(i,1);
		x = y_predict(i,2);
		s = y_predict(i,3);
		DeltaT = T - Tprev;
		asyn_predict(i) = p - DeltaT;
		f_predict = a*(p-(T-s)-baseline) + b*(x-T) ...
			+ alpha*(p-(T-s)-baseline)^3 ...
			+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
			+ gamma*(x-T)^3 ...
			+ baseline;
		g_predict = c*(p-(T-s)-baseline) + d*(x-T) ...
			+ delta*(p-(T-s)-baseline)^2 ...
			+ T;
		h_predict = T;
		y_predict(i+1,1) = f_predict;
		y_predict(i+1,2) = g_predict;
		y_predict(i+1,3) = h_predict;

		% model by hand
		e = y_byhand(i,1);
		x = y_byhand(i,2);
		s = y_byhand(i,3);
		DeltaT = T - Tprev;
		e = e - DeltaT;
		if DeltaT>0
			e_minus(i) = i;
		elseif DeltaT<0
			e_plus(i) = i;
		end
		f_byhand = a*(e-baseline) + b*(x-T) ...
			+ alpha*(e-baseline)^3 ...
			+ beta*(e-baseline)*(x-T)^2 ...
			+ gamma*(x-T)^3 ...
			+ baseline;
		g_byhand = c*(e-baseline) + d*(x-T) ...
			+ delta*(e-baseline)^2 ...
			+ T;
		h_byhand = T;
		asyn_byhand(i) = e;
		y_byhand(i+1,1) = f_byhand;
		y_byhand(i+1,2) = g_byhand;
		y_byhand(i+1,3) = h_byhand;

		% naive model
		e = y_naive(i,1);
		x = y_naive(i,2);
		s = y_naive(i,3);
		DeltaT = T - Tprev;
		f_naive = a*(e-baseline) + b*(x-T) ...
			+ alpha*(e-baseline)^3 ...
			+ beta*(e-baseline)*(x-T)^2 ...
			+ gamma*(x-T)^3 ...
			+ baseline;
		g_naive = c*(e-baseline) + d*(x-T) ...
			+ delta*(e-baseline)^2 ...
			+ T;
		h_naive = T;
		asyn_naive(i) = e;
		y_naive(i+1,1) = f_naive;
		y_naive(i+1,2) = g_naive;
		y_naive(i+1,3) = h_naive;
	end
	
	Tprev = T_all(end-1);
	DeltaT = T_all(end) - Tprev;
	asyn_predict(end) = f_predict - DeltaT;
	asyn_naive(end) = f_naive;
	e = f_byhand - DeltaT;
	asyn_byhand(end) = e;

	if perturb==1
		subplot(4,3,1);
		plot([0 0],y_lims_step_period([1 end]),'k-');
		hold on;
		plot(datos_time_c,T_all,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_step);
		ylim(y_lims_step_period);
		ylabel({'Stimulus';'period T_n (ms)'},'fontsize',fsize,'fontweight','bold','rotation',0,'position',[-13 455]);
		title('Step change T_n','fontsize',fsize,'fontweight','bold');
		text(-8,534,'(a)','fontsize',fsize_text);
		subplot(4,3,4);
		plot([0 0],y_lims_step([1 end]),'k-');
		hold on;
		plot(x_lims_step([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_predict,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'xticklabel',[],'ytick',[0 25 50]);
		xlim(x_lims_step);
		ylim(y_lims_step);
		ylabel({'Model';'response';'e_n (ms)'},'fontsize',fsize,'fontweight','bold','rotation',0,'position',[-13 -10]);
		subplot(4,3,7);
		plot([0 0],y_lims_step([1 end]),'k-');
		hold on;
		plot(x_lims_step([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_byhand,'.-','color',new_cmap(data_color,:),'markersize',msize);
		ptm1 = plot(e_minus-perturb_bip,-10,'o','color',new_cmap(delta_color,:));
		set(ptm1,'markerfacecolor','w','markersize',msize2);
		ptm2 = plot(e_plus-perturb_bip,-10,'o','color',new_cmap(delta_color,:));
		set(ptm2,'markerfacecolor',new_cmap(delta_color,:),'markersize',msize2);
		set(gca,'fontsize',fsize,'xticklabel',[],'ytick',[0 25 50]);
		xlim(x_lims_step);
		ylim(y_lims_step);
		ylabel({'By hand';'e_n (ms)'},'fontsize',fsize,'fontweight','bold','rotation',0,'position',[-13 0]);
		plot(10,20,'o','color',new_cmap(delta_color,:),'markerfacecolor','w','markersize',msize2);
		plot(10,40,'o','color',new_cmap(delta_color,:),'markerfacecolor',new_cmap(delta_color,:),'markersize',msize2);
		text(11,14,'e_n by hand (-)','fontsize',fsize_legend);
		text(11,34,'e_n by hand (+)','fontsize',fsize_legend);
		subplot(4,3,10);
		plot([0 0],y_lims_step([1 end]),'k-');
		hold on;
		plot(x_lims_step([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_naive,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'ytick',[0 25 50]);
		xlim(x_lims_step);
		ylim(y_lims_step);
		xlabel('Step n relative to perturbation');
		ylabel({'Naive';'e_n (ms)'},'fontsize',fsize,'fontweight','bold','rotation',0,'position',[-13 0]);
	elseif perturb==2
		subplot(4,3,2);
		plot([0 0],y_lims_sin_period([1 end]),'k-');
		hold on;
		plot(datos_time_c,T_all,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_sin_rand);
		ylim(y_lims_sin_period);
		title('Sinusoidal T_n','fontsize',fsize,'fontweight','bold');
		text(-9,585,'(b)','fontsize',fsize_text);
		subplot(4,3,5);
		plot([0 0],y_lims_sin([1 end]),'k-');
		hold on;
		plot(x_lims_sin_rand([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_predict,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_sin_rand);
		ylim(y_lims_sin);
		subplot(4,3,8);
		plot([0 0],y_lims_sin([1 end]),'k-');
		hold on;
		plot(x_lims_sin_rand([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_byhand,'.-','color',new_cmap(data_color,:),'markersize',msize);
		ptm = plot(e_minus-perturb_bip,-80,'o','color',new_cmap(delta_color,:));
		set(ptm,'markerfacecolor','w','markersize',msize2);
		ptm = plot(e_plus-perturb_bip,-80,'o','color',new_cmap(delta_color,:));
		set(ptm,'markerfacecolor',new_cmap(delta_color,:),'markersize',msize2);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_sin_rand);
		ylim(y_lims_sin);
		subplot(4,3,11);
		plot([0 0],y_lims_sin([1 end]),'k-');
		hold on;
		plot(x_lims_sin_rand([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_naive,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize);
		xlim(x_lims_sin_rand);
		ylim(y_lims_sin);
		xlabel('Step n relative to perturbation');
	else
		subplot(4,3,3);
		plot([0 0],y_lims_rand_period([1 end]),'k-');
		hold on;
		plot(datos_time_c,T_all,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_sin_rand);
		ylim([450 550]);
		title('Random T_n','fontsize',fsize,'fontweight','bold');
		text(-9,585,'(c)','fontsize',fsize_text);
		subplot(4,3,6);
		plot([0 0],y_lims_rand([1 end]),'k-');
		hold on;
		plot(x_lims_sin_rand([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_predict,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_sin_rand);
		ylim(y_lims_rand);
		subplot(4,3,9);
		plot([0 0],y_lims_rand([1 end]),'k-');
		hold on;
		plot(x_lims_sin_rand([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_byhand,'.-','color',new_cmap(data_color,:),'markersize',msize);
		ptm = plot(e_minus-perturb_bip,-55,'o','color',new_cmap(delta_color,:));
		set(ptm,'markerfacecolor','w','markersize',msize2);
		ptm = plot(e_plus-perturb_bip,-55,'o','color',new_cmap(delta_color,:));
		set(ptm,'markerfacecolor',new_cmap(delta_color,:),'markersize',msize2);
		set(gca,'fontsize',fsize,'xticklabel',[]);
		xlim(x_lims_sin_rand);
		ylim(y_lims_rand);
		subplot(4,3,12);
		plot([0 0],y_lims_rand([1 end]),'k-');
		hold on;
		plot(x_lims_sin_rand([1 end]),[0 0],'k-');
		plot(datos_time_c,asyn_naive,'.-','color',new_cmap(data_color,:),'markersize',msize);
		set(gca,'fontsize',fsize);
		xlim(x_lims_sin_rand);
		ylim(y_lims_rand);
		xlabel('Step n relative to perturbation');
	end

end


filename = 'figure5_benefits';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 6: predicted results (I)


a = 0.981;
b = 0.266;
c = -0.823;
d = 0.0238;
alpha = -2.213e-05;
beta = -7.845e-05;
gamma = 5.338e-05;
delta = 3.348e-3;
T0 = 500;		% ms


lwidth = 1;
msize = 8;
lwidth2 = 2;
msize2 = 3;
fsize = 8;
fsize_legend = 6;
fsize_text = fsize;
fig_size_cm = [20 8];

x_lims = [-5 20];
y_lims_neg = [-30 60];
y_lims_pos = [-80 40];
x_lims_embed = [-90 90];
y_lims_embed = [-90 90];
perturb_bip = 12;
embed_ini = perturb_bip;
embed_fin = perturb_bip + 15;

new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);

figure(6);
clf(6);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);



y_ini = [0 T0 T0];
cmap_color_pos = 1;
cmap_color_neg = 4;
for panel = 1:3
% panel = 2;
	% model parameters
	switch panel
		case 1
			perturb_bip_spatial_pos = perturb_bip + 2;
			perturb_bip_spatial_neg = perturb_bip + 5;
			pert_size_spatial_pos = 0;
			pert_size_spatial_neg = 0;
			pert_size_neg = 4;
			pert_size_pos = 1;
			panel_label = '(a)';
		case 2
			perturb_bip_spatial_pos = perturb_bip + 8;
			perturb_bip_spatial_neg = perturb_bip + 5;
			pert_size_spatial_pos = -40;
			pert_size_spatial_neg = -40;
			pert_size_neg = 4;
			pert_size_pos = 1;
			panel_label = '(b)';
		case 3
			perturb_bip_spatial_pos = perturb_bip + 2;
			perturb_bip_spatial_neg = perturb_bip + 5;
			pert_size_spatial_pos = -40;
			pert_size_spatial_neg = -40;
			panel_label = '(c)';
	end


	% postbaseline set to zero
	perturb_sizes_nobase = [70 60 50 40 30 20 10];
	nbr_perturb_sizes = length(perturb_sizes_nobase);
	% variables
	y = nan(2,nbr_perturb_sizes,nbr_steps,3);
	asyn_nobase = nan(2,nbr_perturb_sizes,nbr_steps);
	% parameter
	T_all = T0*ones(2,nbr_perturb_sizes,nbr_steps);
	T_all(1,:,perturb_bip:end) = T_all(1,:,perturb_bip:end) - repmat(perturb_sizes_nobase,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
	T_all(2,:,perturb_bip:end) = T_all(2,:,perturb_bip:end) + repmat(perturb_sizes_nobase,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
	for pert_sign = 1:2
		for pert_size = 1:nbr_perturb_sizes
			y(pert_sign,pert_size,1,:) = y_ini;
			perturb_size = perturb_sizes_nobase(pert_size);
			baseline = 0;
			for i = 1:nbr_steps-1
				p = y(pert_sign,pert_size,i,1);
				x = y(pert_sign,pert_size,i,2);
				s = y(pert_sign,pert_size,i,3);
				T = T_all(pert_sign,pert_size,i);
				if i==1
					Tprev = T0;
				else
					Tprev = T_all(pert_sign,pert_size,i-1);
				end
				if pert_sign==1
					if i==perturb_bip_spatial_neg
						p = p + pert_size_spatial_neg;
					end
				else
					if i==perturb_bip_spatial_pos
						p = p + pert_size_spatial_pos;
					end
				end
				f = a*(p-(T-s)-baseline) + b*(x-T) ...
					+ alpha*(p-(T-s)-baseline)^3 ...
					+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
					+ gamma*(x-T)^3 ...
					+ baseline;
				g = c*(p-(T-s)-baseline) + d*(x-T) ...
					+ delta*(p-(T-s)-baseline)^2 ...
					+ T;
				h = T;

				DeltaT = T - Tprev;
				asyn_nobase(pert_sign,pert_size,i) = p - DeltaT;

				y(pert_sign,pert_size,i+1,1) = f;
				y(pert_sign,pert_size,i+1,2) = g;
				y(pert_sign,pert_size,i+1,3) = h;
			end
			asyn_nobase(pert_sign,pert_size,end) = p - DeltaT;
		end
	end
	y(y>1000) = NaN;
	asyn_nobase(asyn_nobase>1000) = NaN;
	x_nobase = y(:,:,:,2);
	x_nobase_shift = x_nobase - T_all;



	pert_size = 1;
	subplot(9,3,[panel panel+6]);
	plot([0 0],[y_lims_pos(1) y_lims_neg(2)],'k-');
	hold on;
	plot(x_lims,[0 0],'k-');
	pert_sign = 1;
	legend_text1 = ['-' num2str(perturb_sizes_nobase(pert_size_neg))];
	ptp1 = plot(datos_time,squeeze(asyn_nobase(pert_sign,pert_size_neg,:)),'.-','color',new_cmap(cmap_color_neg,:),'linewidth',lwidth,'markersize',msize);
	pert_sign = 2;
	legend_text2 = ['+' num2str(perturb_sizes_nobase(pert_size_pos))];
	ptp2 = plot(datos_time,squeeze(asyn_nobase(pert_sign,pert_size_pos,:)),'.-','color',new_cmap(cmap_color_pos,:),'linewidth',lwidth,'markersize',msize);
	xlim(x_lims);
	ylim([y_lims_pos(1) y_lims_neg(2)]);
	set(gca,'fontsize',fsize);
	set(gca,'xtick',[x_lims(1):5:x_lims(2)]);
	set(gca,'ytick',[-50:25:50]);
	xlabel('Step n relative to perturbation');
	ylabel('e_n (ms)');
	if panel==1
		ptl = legend([ptp1 ptp2],legend_text1,legend_text2,'location','southeast');
		set(ptl,'fontsize',fsize_legend);
		legend boxoff;
		posf = ds2nfu(gca,[10.2 26 -1.8 -15]);
		pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string','(B)');
		set(pta,'headlength',3,'headwidth',3,'headstyle','plain','linewidth',1,'color',new_cmap(cmap_color_pos,:),'fontsize',0.9*fsize);
		posf = ds2nfu(gca,[3 -25 -0.7 26]);
		pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string','(C)');
		set(pta,'headlength',3,'headwidth',3,'headstyle','plain','linewidth',1,'color',new_cmap(cmap_color_pos,:),'fontsize',0.9*fsize);
		posf = ds2nfu(gca,[6 -25 -1 26]);
		pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string','(B,C)');
		set(pta,'headlength',3,'headwidth',3,'headstyle','plain','linewidth',1,'color',new_cmap(cmap_color_neg,:),'fontsize',0.9*fsize);
	elseif panel==2
		plot(8,-48,'p','color',new_cmap(cmap_color_pos,:),'markerfacecolor',new_cmap(cmap_color_pos,:),'markersize',0.5*msize);
		plot(5,-48,'p','color',new_cmap(cmap_color_neg,:),'markerfacecolor',new_cmap(cmap_color_neg,:),'markersize',0.5*msize);
	elseif panel==3
		plot(2.5,-46,'p','color',new_cmap(cmap_color_pos,:),'markerfacecolor',new_cmap(cmap_color_pos,:),'markersize',0.5*msize);
		plot(5,-47,'p','color',new_cmap(cmap_color_neg,:),'markerfacecolor',new_cmap(cmap_color_neg,:),'markersize',0.5*msize);
	end
	text(-11.5,55,panel_label,'fontsize',fsize_text);


	perturb_sizes_bound_pos_idx = [pert_size_pos];
	perturb_sizes_bound_neg_idx = [pert_size_neg];
	subplot(9,3,[panel+12 panel+24]);
	plot([0 0],y_lims_embed,'k-');
	hold on;
	plot(x_lims_embed,[0 0],'k-');
	pert_sign = 1;
	cmap_color = 4;
	for pert_size_idx = perturb_sizes_bound_neg_idx
		ptp1 = plot(squeeze(asyn_nobase(pert_sign,pert_size_idx,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size_idx,embed_ini:embed_fin)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	end
	pert_sign = 2;
	cmap_color = 1;
	for pert_size_idx = perturb_sizes_bound_pos_idx
		ptp3 = plot(squeeze(asyn_nobase(pert_sign,pert_size_idx,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size_idx,embed_ini:embed_fin)),'.-','color',new_cmap(cmap_color,:),'linewidth',lwidth,'markersize',msize);
	end
	xlim(x_lims_embed);
	ylim(y_lims_embed);
	set(gca,'fontsize',fsize);
	set(gca,'xtick',[-70:35:70]);
	set(gca,'ytick',[-70:35:70]);
	xlabel('e_n (ms)');
	ylabel('x_n - T_{post} (ms)');
	if panel==1
		posf = ds2nfu(gca,[15 8 -5 -15]);
		pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string','(B)');
		set(pta,'headlength',3,'headwidth',3,'headstyle','plain','linewidth',1,'color',new_cmap(cmap_color_pos,:),'fontsize',0.9*fsize);
		posf = ds2nfu(gca,[-5 60 9 12]);
		pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string','(C)');
		set(pta,'headlength',3,'headwidth',3,'headstyle','plain','linewidth',1,'color',new_cmap(cmap_color_pos,:),'fontsize',0.9*fsize);
		posf = ds2nfu(gca,[-9 -21 15 5]);
		pta = annotation('textarrow',[posf(1) posf(1)+posf(3)],[posf(2) posf(2)+posf(4)],'string','(B,C)');
		set(pta,'headlength',3,'headwidth',3,'headstyle','plain','linewidth',1,'color',new_cmap(cmap_color_neg,:),'fontsize',0.9*fsize);
	elseif panel==2
		plot(-40,-14,'p','color',new_cmap(cmap_color_pos,:),'markerfacecolor',new_cmap(cmap_color_pos,:),'markersize',0.5*msize);
		plot(-35,-17,'p','color',new_cmap(cmap_color_neg,:),'markerfacecolor',new_cmap(cmap_color_neg,:),'markersize',0.5*msize);
	elseif panel==3
		plot(-35,69,'p','color',new_cmap(cmap_color_pos,:),'markerfacecolor',new_cmap(cmap_color_pos,:),'markersize',0.5*msize);
		plot(-37,-16,'p','color',new_cmap(cmap_color_neg,:),'markerfacecolor',new_cmap(cmap_color_neg,:),'markersize',0.5*msize);
	end
end


filename = 'figure6_predicted1';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 7: predicted results (II)


a = 0.981;
b = 0.266;
c = -0.823;
d = 0.0238;
alpha = -2.213e-05;
beta = -7.845e-05;
gamma = 5.338e-05;
delta = 3.348e-3;
T0 = 500;		% ms


lwidth = 1;
msize = 15;
lwidth2 = 2;
msize2 = 5;
fsize = 13;
fsize_legend = 10;
fsize_text = fsize;
fig_size_cm = [15 20];

x_lims = [-5 15];
y_lims = [-70 70];
x_lims_embed = [-70 70];
y_lims_embed = [-50 70];
perturb_bip = 12;
embed_ini = perturb_bip;
embed_fin = perturb_bip + 15;

% new_cmap = morgenstemning(9);
% new_cmap = new_cmap(3:end-2,:);
new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);

figure(7);
clf(7);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);



y_ini = [0 T0 T0];
cmap_color = [1 4];
cmap_linestyle = {'-','--'};


% postbaseline set to zero
perturb_sizes_nobase = [60 20];
nbr_perturb_sizes = length(perturb_sizes_nobase);
perturb_signs = [-1 1];
nbr_perturb_signs = length(perturb_signs);
% variables
y = nan(nbr_perturb_signs,nbr_perturb_sizes,nbr_steps,3);
asyn_nobase = nan(nbr_perturb_signs,nbr_perturb_sizes,nbr_steps);
% parameter
T_all = T0*ones(nbr_perturb_signs,nbr_perturb_sizes,nbr_steps);
for pert_sign = 1:nbr_perturb_signs
	perturb_sign = perturb_signs(pert_sign);
	for pert_size = 1:nbr_perturb_sizes
		y(pert_sign,pert_size,1,:) = y_ini;
		perturb_size = perturb_sizes_nobase(pert_size);
		baseline = 0;
		for i = 1:nbr_steps-1
			p = y(pert_sign,pert_size,i,1);
			x = y(pert_sign,pert_size,i,2);
			s = y(pert_sign,pert_size,i,3);
			T = T_all(pert_sign,pert_size,i);
			if i==1
				Tprev = T0;
			else
				Tprev = T_all(pert_sign,pert_size,i-1);
			end
			if i==perturb_bip
				p = p + perturb_sign*perturb_size;
			end
			f = a*(p-(T-s)-baseline) + b*(x-T) ...
				+ alpha*(p-(T-s)-baseline)^3 ...
				+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
				+ gamma*(x-T)^3 ...
				+ baseline;
			g = c*(p-(T-s)-baseline) + d*(x-T) ...
				+ delta*(p-(T-s)-baseline)^2 ...
				+ T;
			h = T;
			
			DeltaT = T - Tprev;
			asyn_nobase(pert_sign,pert_size,i) = p - DeltaT;
			
			y(pert_sign,pert_size,i+1,1) = f;
			y(pert_sign,pert_size,i+1,2) = g;
			y(pert_sign,pert_size,i+1,3) = h;
		end
		asyn_nobase(pert_sign,pert_size,end) = p - DeltaT;
	end
end
y(y>1000) = NaN;
asyn_nobase(asyn_nobase>1000) = NaN;
x_nobase = y(:,:,:,2);
x_nobase_shift = x_nobase - T_all;


subplot(9,1,[1 3]);
plot([0 0],[y_lims_pos(1) y_lims_neg(2)],'k-');
hold on;
plot(x_lims,[0 0],'k-');
for pert_sign = 1:nbr_perturb_signs
	for pert_size = 1:nbr_perturb_sizes
		ptp1 = plot(datos_time,squeeze(asyn_nobase(pert_sign,pert_size,:)),'.','linestyle',cmap_linestyle{pert_sign},'color',new_cmap(cmap_color(pert_size),:),'linewidth',lwidth,'markersize',msize);
	end
end
xlim(x_lims);
ylim(y_lims);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):5:x_lims(2)]);
set(gca,'ytick',[-50:25:50]);
xlabel('Step n relative to perturbation');
ylabel('e_n (ms)');
text(-7.8,60,'(a)','fontsize',fsize_text);


subplot(9,1,[5 9]);
plot([0 0],y_lims_embed,'k-');
hold on;
plot(x_lims_embed,[0 0],'k-');
ptp = [];
for pert_sign = 1:nbr_perturb_signs
	for pert_size = 1:nbr_perturb_sizes
		ptp(pert_sign,pert_size) = plot(squeeze(asyn_nobase(pert_sign,pert_size,embed_ini:embed_fin)),squeeze(x_nobase_shift(pert_sign,pert_size,embed_ini:embed_fin)),'.','linestyle',cmap_linestyle{pert_sign},'color',new_cmap(cmap_color(pert_size),:),'linewidth',lwidth,'markersize',msize);
	end
end
xlim(x_lims_embed);
ylim(y_lims_embed);
set(gca,'fontsize',fsize);
set(gca,'xtick',[-60:30:60]);
set(gca,'ytick',[-60:30:60]);
xlabel('e_n (ms)');
ylabel('x_n - T_{post} (ms)');
legend([ptp(1,1) ptp(2,1) ptp(1,2) ptp(2,2)],['-' num2str(perturb_sizes_nobase(1)) ' ms'],['+' num2str(perturb_sizes_nobase(1)) ' ms'],['-' num2str(perturb_sizes_nobase(2)) ' ms'],['+' num2str(perturb_sizes_nobase(2)) ' ms'],'location','northeast');
legend boxoff;
text(-90,65,'(b)','fontsize',fsize_text);



filename = 'figure7_predicted2';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);




%% Figure 8: return points


%load data
data_filename = 'data_en_nopostbase_zn.dat';
datos = load(data_filename);

% time series
datos_asyn_b = [];
datos_time_b = datos(:,1);
datos_asyn_b(1,:,:,1) = datos(:,[2:2:10]);		% negative steps en
datos_asyn_b(2,:,:,1) = fliplr(datos(:,[12:2:end]));	% positive steps en
datos_asyn_b(1,:,:,2) = datos(:,[3:2:11]);		% negative steps zn
datos_asyn_b(2,:,:,2) = fliplr(datos(:,[13:2:end]));	% positive steps zn


lwidth = 1;
msize = 10;%8;
lwidth2 = 1;
%msize2 = 3;
fsize = 10;
fsize_legend = 10;
fsize_text = 10;
fig_size_cm = [25 10];



x_lims = [-2 16];
y_lims_neg = [-10 50];
y_lims_pos = [-50 20];
y_lims_neg_z = [-20 40];
y_lims_pos_z = [-40 30];

new_cmap = morgenstemning(9);
new_cmap = new_cmap(3:end-2,:);


a1=-1.1125;
b1=10.275;
c1=-32.5875;
d1=29.325;
e1=34.822;

a2 = 0.597917;
b2 = -8.15083;
c2 = 33.4671;
d2 = -33.4142;
e2 = -36.7686;

a3 = -1.84696;
b3 = 14.0196;
c3 = -27.1435;
d3 = -8.10208;
e3 = 28.973;
       
a4 = 1.85421;
b4 = -16.1504;
c4 = 36.8573;
d4 = 1.06792;
e4 = -31.129;

 

figure(12);
clf(12);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

pert_sign = 1;
ptp = [];
subplot(2,2,1);
plot([0 0],y_lims_neg,'k-');
hold on;
plot(x_lims,[0 0],'k-');
pert_size =3;
ptp(1) = plot(datos_time_b,squeeze(datos_asyn_b(pert_sign,:,pert_size,1)),'.--','color',new_cmap(4,:),'linewidth',lwidth,'markersize',msize);

x1 = linspace(-0.2,4.2);
p = [a1 b1 c1 d1 e1];
y1=polyval(p,x1);
ptp(2)=plot(x1,y1,'-','color',new_cmap(1,:),'linewidth',lwidth2);
set(gca,'FontSize',fsize);
text(-5,50,'(a)','FontSize',fsize);

xlim(x_lims);
ylim(y_lims_neg);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):2:x_lims(2)]);
set(gca,'ytick',[y_lims_neg(1):10:y_lims_neg(2)]+10);
ylabel('e_n (ms)');
ptl = legend(ptp,'-30','fit','location','northeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;


pert_sign = 2;
ptp = [];
subplot(2,2,2);
plot([0 0],y_lims_pos,'k-');
hold on;
plot(x_lims,[0 0],'k-');
pert_size =3;
ptp(1)= plot(datos_time_b,squeeze(datos_asyn_b(pert_sign,:,pert_size,1)),'.--','color',new_cmap(4,:),'linewidth',lwidth,'markersize',msize);

x1 = linspace(-0.2,4.2);
p = [a2 b2 c2 d2 e2];
y1=polyval(p,x1);
ptp(2)=plot(x1,y1,'-','color',new_cmap(1,:),'linewidth',lwidth2);
set(gca,'FontSize',fsize);
text(-5,20,'(b)','FontSize',fsize);

xlim(x_lims);
ylim(y_lims_pos);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):2:x_lims(2)]);
set(gca,'ytick',[y_lims_pos(1):10:y_lims_pos(2)]+10);
ylabel('e_n (ms)');
ptl = legend(ptp,'+30','fit','location','southeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;

pert_sign = 1;
ptp = [];
subplot(2,2,3);
plot([0 0],y_lims_neg_z,'k-');
hold on;
plot(x_lims,[0 0],'k-');
pert_size =3;
ptp(1)= plot(datos_time_b,squeeze(datos_asyn_b(pert_sign,:,pert_size,2)),'.--','color',new_cmap(4,:),'linewidth',lwidth,'markersize',msize);

x1 = linspace(-0.2,4.2);
p = [a3 b3 c3 d3 e3];
y1=polyval(p,x1);
ptp(2)=plot(x1,y1,'-','color',new_cmap(1,:),'linewidth',lwidth2);

xlim(x_lims);
ylim(y_lims_neg_z);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):2:x_lims(2)]);
set(gca,'ytick',[y_lims_neg_z(1):10:y_lims_neg_z(2)]+10);
xlabel('Step n relative to perturbation');
ylabel('z_n = e_n - e_n_-_1 (ms)');
ptl = legend(ptp,'-30','fit','location','northeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;

pert_sign = 2;
ptp = [];
subplot(2,2,4);
plot([0 0],y_lims_pos_z,'k-');
hold on;
plot(x_lims,[0 0],'k-');
pert_size =3;
ptp(1) = plot(datos_time_b,squeeze(datos_asyn_b(pert_sign,:,pert_size,2)),'.--','color',new_cmap(4,:),'linewidth',lwidth,'markersize',msize);

x1 = linspace(-0.2,4.2);
p = [a4 b4 c4 d4 e4];
y1=polyval(p,x1);
ptp(2)=plot(x1,y1,'-','color',new_cmap(1,:),'linewidth',lwidth2);

xlim(x_lims);
ylim(y_lims_pos_z);
set(gca,'fontsize',fsize);
set(gca,'xtick',[x_lims(1):2:x_lims(2)]);
set(gca,'ytick',[y_lims_pos_z(1):10:y_lims_pos_z(2)]+10);
xlabel('Step n relative to perturbation');
ylabel('z_n = e_n - e_n_-_1  (ms)');
ptl = legend(ptp,'+30','fit','location','southeast');
set(ptl,'fontsize',fsize_legend);
legend boxoff;


headWidth = 6; 
headLength = 8;
headWidth2 = 4; 
headLength2 = 6;

% Create arrow from emax to z(emax)
annotation(figure(12),'arrow',[0.17875 0.17875],...
   [0.882517482517483 0.331468531468532],'LineWidth',1, 'Linestyle','--','HeadLength',headLength,'HeadWidth',headWidth);

      
% Create textarrow for indicate emax
annotation(figure(12),'textarrow',[0.215625 0.181875],...
    [0.888 0.888],'LineWidth',1,'headStyle','plain','HeadLength',headLength2,'HeadWidth',headWidth2);

% Create textbox for emax
annotation(figure(12),'textbox',...
    [0.215 0.86 0.034375 0.0521468531468531],...
    'String',{'e_{ret}'},...
    'LineStyle','none',...
    'FontSize',9,...
    'FitBoxToText','off');

% Create textarrow for indicate z(emax)
annotation(figure(12),'textarrow',[0.21625 0.181875],...
    [0.331468531468532 0.331468531468532],'LineWidth',1,'headStyle','plain','HeadLength',headLength2,'HeadWidth',headWidth2);
    
% Create textbox for z(emax)
annotation(figure(12),'textbox',...
    [0.215 0.27 0.11281249682419 0.0843146853146853],...
    'String',{'z(e_{ret})'},...
    'LineStyle','none',...
    'FontSize',9,...
    'FitBoxToText','off');
    
% Create arrow from zmax to e(zmax)
annotation(figure(12),'arrow',[0.6531274 0.6531274],...
    [0.410188811188811 0.746],'LineWidth',1,'Linestyle','--','HeadLength',headLength,'HeadWidth',headWidth);

% Create textarrow for e(zmax)
annotation(figure(12),'textarrow',[0.69253125 0.65815625],...
    [0.746 0.746],'LineWidth',1,'headStyle','plain','HeadLength',headLength2,'HeadWidth',headWidth2);

% Create textarrow for zmax
annotation(figure(12),'textarrow',[0.69253125 0.65815625],...
    [0.412280701754386 0.412280701754386],'LineWidth',1,'headStyle','plain','HeadLength',headLength2,'HeadWidth',headWidth2);


% Create textbox for e(zmax)
annotation(figure(12),'textbox',...
    [0.69 0.71 0.0742187480209396 0.0628654956381921],...
    'String',{'e(z_{ret})'},...
    'LineStyle','none',...
    'FontSize',9);
   
   
% Create textbox for zmax
annotation(figure(12),'textbox',...
    [0.69 0.37 0.0515624986961485 0.0647894736842106],...
    'String',{'z_{ret}'},...
    'LineStyle','none',...
    'FontSize',9,...
    'FitBoxToText','off');



filename = 'figure8_return';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);




%% Figure 9: all perturbation sizes


% model parameters
a = 0.981;
b = 0.266;
c = -0.823;
d = 0.0238;
alpha = -2.213e-05;
beta = -7.845e-05;
gamma = 5.338e-05;
delta = 3.348e-3;
T0 = 500;		% ms

postbl = squeeze(mean(datos_asyn(:,postbl_ini:end,:),2));
exp_postbaseline = postbl';

nbr_steps = length(datos_time);
perturb_bip = 12;
y_ini = [0 T0 T0];


% postbaseline set to experimental values
perturb_sizes = [50 40 30 20 10];
nbr_perturb_sizes = length(perturb_sizes);
% variables
y = nan(2,nbr_perturb_sizes,nbr_steps,3);
asyn = nan(2,nbr_perturb_sizes,nbr_steps);
% parameter
T_all = T0*ones(2,nbr_perturb_sizes,nbr_steps);
T_all(1,:,perturb_bip:end) = T_all(1,:,perturb_bip:end) - repmat(perturb_sizes,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
T_all(2,:,perturb_bip:end) = T_all(2,:,perturb_bip:end) + repmat(perturb_sizes,[1 1 size(T_all(1,:,perturb_bip:end),3)]);
for pert_sign = 1:2
	for pert_size = 1:nbr_perturb_sizes
		y(pert_sign,pert_size,1,:) = y_ini;
		perturb_size = perturb_sizes(pert_size);
		baseline = 0;
		for i = 1:nbr_steps-1
			p = y(pert_sign,pert_size,i,1);
			x = y(pert_sign,pert_size,i,2);
			s = y(pert_sign,pert_size,i,3);
			T = T_all(pert_sign,pert_size,i);
			if i==1
				Tprev = T0;
			else
				Tprev = T_all(pert_sign,pert_size,i-1);
			end
			if i==perturb_bip
				baseline = exp_postbaseline(pert_size,pert_sign);
			end
			f = a*(p-(T-s)-baseline) + b*(x-T) ...
				+ alpha*(p-(T-s)-baseline)^3 ...
				+ beta*(p-(T-s)-baseline)*(x-T)^2 ...
				+ gamma*(x-T)^3 ...
				+ baseline;
			g = c*(p-(T-s)-baseline) + d*(x-T) ...
				+ delta*(p-(T-s)-baseline)^2 ...
				+ T;
			h = T;
			
			DeltaT = T - Tprev;
			asyn(pert_sign,pert_size,i) = p - DeltaT;
			
			y(pert_sign,pert_size,i+1,1) = f;
			y(pert_sign,pert_size,i+1,2) = g;
			y(pert_sign,pert_size,i+1,3) = h;
		end
		asyn(pert_sign,pert_size,end) = p - DeltaT;
	end
end
y(y>1000) = NaN;
asyn(asyn>1000) = NaN;
x = y(:,:,:,2);
x_shift = x - T_all;





lwidth = 1;
msize2 = 11;
lwidth2 = 2;
msize1 = 4;
fsize = 11;
fsize_legend = 10;
fsize_text = fsize;
panel_labels = {'(a)','(b)','(c)','(d)','(e)'};
fig_size_cm = [30 12];

embed_ini = 11;
embed_fin = 20;
x_lims = [-5 20];
y_lims = [-60 60];

new_cmap = morgenstemning(10);
new_cmap = new_cmap([1 5 6 7 8],:);

figure(9);
clf(9);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

for pert_size = 1:nbr_perturb_sizes
% 	subplot(1,nbr_perturb_sizes,pert_size);
	subplot(2,3,pert_size);
	pert_sign = 1;
	cmap = 4;
	plot([0 0],y_lims,'k-');
	hold on;
	plot(x_lims,[0 0],'k-');
	ptp_exp_neg = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'o-','color',new_cmap(cmap,:),'linewidth',lwidth,'markersize',msize1);
	ptp_leg_neg = plot(datos_time,squeeze(asyn(pert_sign,pert_size,:)),'-','color',new_cmap(cmap,:),'linewidth',lwidth,'markersize',msize2);
	ptp_fit_neg = plot(datos_time,squeeze(asyn(pert_sign,pert_size,:)),'.-','color',new_cmap(cmap,:),'linewidth',lwidth,'markersize',msize2);
	pert_sign = 2;
	cmap = 1;
	ptp_exp_pos = plot(datos_time,squeeze(datos_asyn(pert_sign,:,pert_size)),'o-','color',new_cmap(cmap,:),'linewidth',lwidth,'markersize',msize1);
	ptp_leg_pos = plot(datos_time,squeeze(asyn(pert_sign,pert_size,:)),'-','color',new_cmap(cmap,:),'linewidth',lwidth,'markersize',msize2);
	ptp_fit_pos = plot(datos_time,squeeze(asyn(pert_sign,pert_size,:)),'.-','color',new_cmap(cmap,:),'linewidth',lwidth,'markersize',msize2);

	xlim(x_lims);
	ylim(y_lims);
	set(gca,'fontsize',fsize);
	set(gca,'xtick',[-5:5:20]);
	set(gca,'ytick',[-50:25:50]);
	xlabel('Step n relative to perturbation');
	ylabel('e_n (ms)');
	ptl = legend([ptp_leg_neg ptp_leg_pos],['-' num2str(perturb_sizes(pert_size))],['+' num2str(perturb_sizes(pert_size))],'location','southeast');
	set(ptl,'fontsize',fsize_legend);
	legend boxoff;
	if pert_size==1
		text(15,45,'o  exp','fontsize',fsize_legend);
		text(15,35,'\bullet  fit','fontsize',fsize_legend);
	end
	text(-10.5,55,panel_labels{pert_size},'fontsize',fsize_text);
end



filename = 'figure9_allfits';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 10: parameter distribution


input_file = 'evols_GBL.dat';
input_file_fid = fopen(input_file);
evols = textscan(input_file_fid,'%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f','HeaderLines',1);

param_names_full = {'a','b','c','d','alfa1','beta1','gamma1','delta1','epsilon1','dseta1','eta1','alfa2','beta2','gamma2','delta2','epsilon2','dseta2','eta2'};
param_names = {'a','b','c','d','delta1','dseta1','eta1','alfa2'};
param_names_latex = {'a','b','c','d','\alpha','\beta','\gamma','\delta'};	% as in the manuscript
param_numbers = find(sum(cell2mat(cellfun(@(x) strcmp(x,param_names_full)',param_names,'UniformOutput',false))'));
param_orders = {'linear','linear','linear','linear','cubic','cubic','cubic','quadratic'};
param_units = {'nondim','nondim','nondim','nondim','ms^{-2}','ms^{-2}','ms^{-2}','ms^{-1}'};
param_units_prefix = {'','','','','10^{-5} ','10^{-5} ','10^{-5} ','10^{-3} '};
param_units_prefix_factor = [1 1 1 1 1e-5 1e-5 1e-5 1e-3];


n_data = size(evols{1},1);
nro_loop = evols{1};
params = [evols{param_numbers+2}];
fitness = evols{2};

[nro_loops,nro_params] = size(params);
[fitness_sort,nro_loop_sort] = sort(fitness,'descend');
params_sort = params(nro_loop_sort,:);


% calcula todos los histogramas individuales
n_bins = 20;
param_hist = zeros(n_bins,nro_params);
param_bin = zeros(n_bins,nro_params);
for i=1:nro_params
	[param_hist(:,i),param_bin(:,i)] = hist(params(:,i),n_bins);
end
hist_lim_x = [[-1 1];...
			[-1 1];...
			[-1 1];...
			[-1 1];...
			[-10 10];...
			[-10 10];...
			[-10 10];...
			[-10 10]];


fsize = 8;
lwidth = 3;
msize = 15;
lwidth2 = 2;
msize2 = 5;
fig_size_cm = [20 8];

new_cmap = morgenstemning(9);
new_cmap = new_cmap(3,:);

figure(8);
clf(8);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

for param = 1:nro_params
	subplot(2,4,param);
	set(gca,'fontsize',fsize);
	ptb = bar(param_bin(:,param)/param_units_prefix_factor(param),param_hist(:,param));
	hold on;
	set(ptb,'facecolor',new_cmap(1,:));
	xlim(hist_lim_x(param,:));
	set(gca,'xtick',[hist_lim_x(param,1):hist_lim_x(param,2)/2:hist_lim_x(param,2)]);
	xlabel([param_names_latex{param} ' (' param_units_prefix{param} param_units{param} ')']);
end


filename = 'figure10_params';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 11: param distribution (linear vs linear)


n_bins_x = 20;
n_bins_y = 20;
n_bins_z = 20;


fsize = 12;
fsize_xlabel = 12;
msize = 20;
lwidth = 3;
fig_size_cm = [15 15];

new_cmap = morgenstemning(9);
new_cmap = new_cmap(3,:);

figure(9);
clf(9);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);
hist2D_colormap = [[1 1 1]; [0 0 1]];
colormap(hist2D_colormap);


for param1_idx = 1:3
	param1 = param1_idx+1;
	data_x = params(:,param1);
	min_x = min(data_x);
	max_x = max(data_x);
	range_x = max_x - min_x;
	bin_x = range_x/n_bins_x;
	x_edges = [min_x-bin_x:bin_x:max_x+bin_x];
	param1_bin_x = x_edges;
	param1_hist = histc(params(:,param1),param1_bin_x);

	subplot(7,7,[(2*param1_idx-1) 2*param1_idx]);
	% distribución individual de los lineales
	ptb1 = bar(x_edges,param1_hist,'histc');
	set(ptb1,'facecolor',new_cmap(1,:));
	set(gca,'xticklabel',[],'ytick',[0:25:50],'fontsize',fsize);
	if param1_idx~=1
		set(gca,'yticklabel',[]);
	end
	xlim(x_edges([1 end]));
	ylim([0 50]);
    
    

	for param2_idx = 1:3
		param2 = param2_idx;
		data_y = params(:,param2);
		min_y = min(data_y);
		max_y = max(data_y);
		range_y = max_y - min_y;
		bin_y = range_y/n_bins_y;
		y_edges = [min_y-bin_y:bin_y:max_y+bin_y];
		param2_bin_x = y_edges;
		param2_hist = histc(params(:,param2),param2_bin_x);

		histmat = hist2(data_x,data_y,x_edges,y_edges);
		
		[aux_x,aux_y] = find(histmat>0);
		histmatBN = double(histmat>0);


		if (param1_idx==1 && param2_idx==3) || (param1_idx==1 && param2_idx==2) || (param1_idx==2 && param2_idx==3)
		else
			subplot(7,7,[(2*param1_idx + (2*param2_idx-1)*7 - 1) (2*param1_idx+7 + (2*param2_idx-1)*7)]);
			% distribución conjunta
			histmatBN(histmatBN==0) = NaN;
			ptp = pcolor(x_edges,y_edges,histmatBN');
			if param1_idx==param2_idx
				xlabel([param_names{param1} ' (' param_units{param1} ')'],'fontsize',fsize_xlabel);
			else
				set(gca,'xtick',[],'xticklabel',[]);
			end
			if param1_idx==param2_idx
				ylabel([param_names{param2} ' (' param_units{param2} ')'],'fontsize',fsize_xlabel);
			else
				set(gca,'ytick',[],'yticklabel',[]);
			end
			set(gca,'fontsize',fsize);
		end

		subplot(7,7,[7*(2*param2_idx) 7*(2*param2_idx)+7]);
		% distribución individual de los nolineales
		ptb3 = barh(y_edges,param2_hist,'histc');
		set(ptb3,'facecolor',new_cmap(1,:));
		if param2_idx~=3
			set(gca,'xticklabel',[]);
		end
		set(gca,'yticklabel',[],'xtick',[0:25:50],'fontsize',fsize);
		xlim([0 50]);
		ylim(y_edges([1 end]));
		
		
	end
end

colormap(new_cmap);

filename = 'figure11_linlin';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);




%% Figure 12: param distribution (nonlinear vs nonlinear)


n_bins_x = 20;
n_bins_y = 20;
n_bins_z = 20;


fsize = 10;
fsize_xlabel = 10;
msize = 20;
lwidth = 3;
fig_size_cm = [15 15];

new_cmap = morgenstemning(9);
new_cmap = new_cmap(3,:);

figure(10);
clf(10);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);
hist2D_colormap = [[1 1 1]; [0 0 1]];
colormap(hist2D_colormap);


for param1_idx = 1:3
	param1 = param1_idx+5;
	data_x = params(:,param1);
	min_x = min(data_x);
	max_x = max(data_x);
	range_x = max_x - min_x;
	bin_x = range_x/n_bins_x;
	x_edges = [min_x-bin_x:bin_x:max_x+bin_x];
	param1_bin_x = x_edges;
	param1_hist = histc(params(:,param1),param1_bin_x);

	subplot(7,7,[(2*param1_idx-1) 2*param1_idx]);
	% distribución individual de los lineales
	ptb1 = bar(x_edges/param_units_prefix_factor(param1),param1_hist,'histc');
	set(ptb1,'facecolor',new_cmap(1,:));
	set(gca,'xticklabel',[],'ytick',[0:40:80],'fontsize',fsize);
	if param1_idx~=1
		set(gca,'yticklabel',[]);
	end
	xlim(x_edges([1 end])/param_units_prefix_factor(param1));
	ylim([0 80]);
       
   
    
	for param2_idx = 1:3
		param2 = param2_idx+4;
		data_y = params(:,param2);
		min_y = min(data_y);
		max_y = max(data_y);
		range_y = max_y - min_y;
		bin_y = range_y/n_bins_y;
		y_edges = [min_y-bin_y:bin_y:max_y+bin_y];
		param2_bin_x = y_edges;
		param2_hist = histc(params(:,param2),param2_bin_x);

		histmat = hist2(data_x,data_y,x_edges,y_edges);
		
		[aux_x,aux_y] = find(histmat>0);
		histmatBN = double(histmat>0);


		if (param1_idx==1 && param2_idx==3) || (param1_idx==1 && param2_idx==2) || (param1_idx==2 && param2_idx==3)
		else
			subplot(7,7,[(2*param1_idx + (2*param2_idx-1)*7 - 1) (2*param1_idx+7 + (2*param2_idx-1)*7)]);
			% distribución conjunta
			histmatBN(histmatBN==0) = NaN;
			ptp = pcolor(x_edges/param_units_prefix_factor(param1),y_edges/param_units_prefix_factor(param2),histmatBN');
			if param1_idx==param2_idx
				xlabel([param_names_latex{param1} ' (' param_units_prefix{param1} ' ' param_units{param1} ')  '],'fontsize',fsize_xlabel);
			else
				set(gca,'xtick',[],'xticklabel',[]);
			end
			if param1_idx==param2_idx
				ylabel([param_names_latex{param2} ' (' param_units_prefix{param2} ' ' param_units{param2} ')  '],'fontsize',fsize_xlabel);
			else
				set(gca,'ytick',[],'yticklabel',[]);
			end
			set(gca,'fontsize',fsize);
		end

		subplot(7,7,[7*(2*param2_idx) 7*(2*param2_idx)+7]);
		% distribución individual de los nolineales
		ptb3 = barh(y_edges/param_units_prefix_factor(param2),param2_hist,'histc');
		set(ptb3,'facecolor',new_cmap(1,:));
		set(gca,'xtick',[0:40:80],'yticklabel',[],'fontsize',fsize);
		if param2_idx~=3
			set(gca,'xticklabel',[]);
		end
		xlim([0 80]);
		ylim(y_edges([1 end])/param_units_prefix_factor(param2));
		
		
	end
end

colormap(new_cmap);


filename = 'figure12_nolinnolin';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);



%% Figure 13: param distribution (linear vs nonlinear)


n_bins_x = 20;
n_bins_y = 20;
n_bins_z = 20;


fsize = 8;
fsize_xlabel = 8;
msize = 20;
lwidth = 3;
fig_size_cm = [15 15];

new_cmap = morgenstemning(9);
new_cmap = new_cmap(3,:);

figure(11);
clf(11);
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);
hist2D_colormap = [[1 1 1]; [0 0 1]];
colormap(hist2D_colormap);


for param1_idx = 1:4
	param1 = param1_idx;
	data_x = params(:,param1);
	min_x = min(data_x);
	max_x = max(data_x);
	range_x = max_x - min_x;
	bin_x = range_x/n_bins_x;
	x_edges = [min_x-bin_x:bin_x:max_x+bin_x];
	param1_bin_x = x_edges;
	param1_hist = histc(params(:,param1),param1_bin_x);

	subplot(9,9,[(2*param1_idx-1) 2*param1_idx]);
	% distribución individual de los lineales
	ptb1 = bar(x_edges,param1_hist,'histc');
	set(ptb1,'facecolor',new_cmap(1,:));
	set(gca,'xticklabel',[],'ytick',[0:25:50],'fontsize',fsize);
	if param1_idx~=1
		set(gca,'yticklabel',[]);
	end
	xlim(x_edges([1 end]));
	ylim([0 50]);

	for param2_idx = 1:4
		param2 = param2_idx + 4;
		data_y = params(:,param2);
		min_y = min(data_y);
		max_y = max(data_y);
		range_y = max_y - min_y;
		bin_y = range_y/n_bins_y;
		y_edges = [min_y-bin_y:bin_y:max_y+bin_y];
		param2_bin_x = y_edges;
		param2_hist = histc(params(:,param2),param2_bin_x);

		histmat = hist2(data_x,data_y,x_edges,y_edges);
		
		[aux_x,aux_y] = find(histmat>0);
		histmatBN = double(histmat>0);


		subplot(9,9,[(2*param1_idx + (2*param2_idx-1)*9 - 1) (2*param1_idx+9 + (2*param2_idx-1)*9)]);
		% distribución conjunta
		histmatBN(histmatBN==0) = NaN;
		ptp = pcolor(x_edges,y_edges/param_units_prefix_factor(param2),histmatBN');
		if param2_idx==4
			xlabel([param_names{param1} ' (' param_units{param1} ')'],'fontsize',fsize_xlabel);
		else
			set(gca,'xtick',[],'xticklabel',[]);
		end
		if param1_idx==1
			ylabel([param_names_latex{param2} ' (' param_units_prefix{param2} ' ' param_units{param2} ')'],'fontsize',fsize_xlabel);
		else
			set(gca,'ytick',[],'yticklabel',[]);
		end
		set(gca,'fontsize',fsize);

		subplot(9,9,[9*(2*param2_idx) 9*(2*param2_idx)+9]);
		% distribución individual de los nolineales
		ptb3 = barh(y_edges/param_units_prefix_factor(param2),param2_hist,'histc');
		set(ptb3,'facecolor',new_cmap(1,:));
		set(gca,'xtick',[0:40:80],'yticklabel',[],'fontsize',fsize);
		if param2_idx~=4
			set(gca,'xticklabel',[]);
		end
		xlim([0 80]);
		ylim(y_edges([1 end])/param_units_prefix_factor(param2));
		
		
	end
end

colormap(new_cmap);

filename = 'figure13_linnolin';
print('-dpdf',[filename '.pdf']);
print('-dpng','-r300',[filename '.png']);




%%

