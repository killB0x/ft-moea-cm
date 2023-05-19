clc, clear all, close all

check_folders = {'data_for_figures/fig9/parametric_spurious_variables/monopropellant_propulsion_system_noise_0per_MCS1TS1Acc1_sfv6';...
    'data_for_figures/fig9/parametric_varying_initial_parent_FTs/parent_disjunctive_normal_form/monopropellant_propulsion_system_noise_0per_MCS1TS1Acc1_sfv6';...
    'data_for_figures/fig9/parametric_varying_initial_parent_FTs/oversized_fault_tree/monopropellant_propulsion_system_noise_0per_MCS1TS1Acc1_sfv6'};

saving_fig = 'figures_in_paper';
save_plot = 1;
leg = {'Setup A','Setup B','Setup C'};

compare = {};
output_results = {};
compare_all = {};
for i = 1 : length(check_folders)
    output_results = [output_results;get_all_results_output(get_dir_files_str(check_folders{i}))];
    compare = [compare;get_all_results(get_dir_files_str(check_folders{i}))];
    compare_all = [compare_all;get_all_all_results(get_dir_files_str(check_folders{i}))];
end

%%
fig_position = [ -281   177   210   193];
axis_position = [0.3019    0.2812    0.6462    0.6786];
text_size= 16;
plot_setup.font_name  = 'Arial Unicode MS';
linewidth = 2;
trial = 1;

figure('color',[1 1 1],'position',fig_position), hold on
for i = 1 : length(check_folders)
    plot(1-compare{i}(trial).acc_data,'marker','.','linestyle','-','linewidth',linewidth)
end
leg1 = legend(leg);
set(leg1,'location','best')
xlabel('Generation')
ylabel('\phi_d')
box on
set(gca,'position',axis_position,'FontName',plot_setup.font_name,'FontSize',text_size,'XTick',[0:10:60])%,'YTick',[0.6:0.1:1])
xlim([0 60])
if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',saving_fig,'/acc_data_vs_generations']),'-dpdf','-r0')
    savefig(gcf,strcat(['',saving_fig,'/acc_data_vs_generations']))
    close(gcf)
end

figure('color',[1 1 1],'position',fig_position), hold on
for i = 1 : length(check_folders)
    plot(1-compare{i}(trial).acc_mcs,'marker','.','linestyle','-','linewidth',linewidth)
end
leg1 = legend(leg);
set(leg1,'location','best')
xlabel('Generation')
ylabel('\phi_c')
box on
set(gca,'position',axis_position,'FontName',plot_setup.font_name,'FontSize',text_size,'XTick',[0:10:60])%,'YTick',[0:0.2:1])
xlim([0 60])
if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',saving_fig,'/acc_mcs_vs_generations']),'-dpdf','-r0')
    savefig(gcf,strcat(['',saving_fig,'/acc_mcs_vs_generations']))
    close(gcf)
end


figure('color',[1 1 1],'position',fig_position), hold on
for i = 1 : length(check_folders)
    plot(compare{i}(trial).tree_size,'marker','.','linestyle','-','linewidth',linewidth)
end
plot([0 max(double(get(gca,'XTick')))+1],[23 23],'r-.')
leg1 = legend(leg);
set(leg1,'location','best')
xlabel('Generation')
ylabel('\phi_s')
box on
set(gca,'position',axis_position,'FontName',plot_setup.font_name,'FontSize',text_size,'XTick',[0:10:60],'YTick',[0:10:50])
xlim([0 60])
ylim([0 50])
if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',saving_fig,'/tree_size_vs_generations']),'-dpdf','-r0')
    savefig(gcf,strcat(['',saving_fig,'/tree_size_vs_generations']))
    close(gcf)
end



figure('color',[1 1 1],'position',fig_position), hold on
for i = 1 : length(check_folders)
    plot(compare{i}(trial).time/60,'marker','.','linestyle','-','linewidth',linewidth)
end
leg1 = legend(leg);
set(leg1,'location','best')
xlabel('Generation')
ylabel('t (min)')
box on
set(gca,'position',axis_position,'FontName',plot_setup.font_name,'FontSize',text_size,'XTick',[0:10:60])
xlim([0 60])
if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',saving_fig,'/time_size_vs_generations']),'-dpdf','-r0')
    savefig(gcf,strcat(['',saving_fig,'/time_size_vs_generations']))
    close(gcf)
end


figure('color',[1 1 1],'position',fig_position), hold on
for i = 1 : length(check_folders)
    plot(diff(compare{i}(trial).time/60),'marker','.','linestyle','-','linewidth',linewidth)
end
leg1 = legend(leg);
set(leg1,'location','best')
xlabel('(To) Generation')
ylabel('\Delta_t (min)')
box on
set(gca,'position',axis_position,'FontName',plot_setup.font_name,'FontSize',text_size,'XTick',[0:10:60])
xlim([0 60])
if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',saving_fig,'/diff_time_size_vs_generations']),'-dpdf','-r0')
    savefig(gcf,strcat(['',saving_fig,'/diff_time_size_vs_generations']))
    close(gcf)
end

%% 
%createfigure(compare_all{1}.tree_size(:,[1,2,3,10,20]), compare_all{1}.acc_data(:,[1,2,3,10,20]), compare_all{1}.acc_mcs(:,[1,2,3,10,20]))
%createfigure(compare_all{2}.tree_size(:,[1,2,3,10,20]), compare_all{2}.acc_data(:,[1,2,3,10,20]), compare_all{2}.acc_mcs(:,[1,2,3,10,20]))
%exportgraphics(gcf,'d_s.png','Resolution',300)
% figure('color',[1 1 1]),hold on
% for i = 1 : size(compare_all{2}.tree_size,2)
%     plot3(compare_all{2}.tree_size(:,i),compare_all{2}.acc_data(:,i),compare_all{2}.acc_mcs(:,i),'*','displayname',strcat(['Gen. ',num2str(i),'']))
% end
%%
function all_files= get_dir_files_str(folder_results)
%% Get the parent directories:
files = dir(folder_results);
files = {files.name}';
a = [];
for i = 1:length(files)
    if files{i}(1) ~= 'c'
        a = [a;i];
    end
end
files(a) = [];
g = [];
for i = 1 : length(files)
    a = strfind(files{i},'_');
    g = [g;str2num(files{i}(a(end)+1:end))];
end
[~,x] = sort(g,'ascend');
s = g(x);
files = files(x);
%% Working directories:
all_files = [];
for i = 1 : length(files)
    dir_path = strcat(['',folder_results,'/',files{i},'/run_0']);
    subfiles = dir(dir_path);
    subfiles = {subfiles.name}';
    a = [];
    for j = 1:length(subfiles)
        if subfiles{j}(1) ~= 'g'
            a = [a;j];
        end
    end
    subfiles(a) = [];
    g = [];
    for j = 1 : length(subfiles)
        a = strfind(subfiles{j},'_');
        g = [g;str2num(subfiles{j}(a(1)+1:end-4))];
    end
    [~,x] = sort(g,'ascend');
    subfiles = subfiles(x);
    %subfiles
    all_files(i).pop_size = s(i);
    all_files(i).path = dir_path;
    all_files(i).files = subfiles;
end
end
% Get results:
function r = get_all_results(c)
for j = 1 : length(c)
    a = c(j);
    pop_size = a.pop_size;
    acc_data = [];
    tree_size =[];
    time = [];
    acc_mcs=[];
    str_fts = {};
    for i = 1 : length(a.files)
        b = load(strcat(['',a.path,'/',a.files{i},'']));
        acc_data = [acc_data;b.parameters_trimmed(end,3)];
        time = [time;b.time];
        tree_size =[tree_size;b.parameters_trimmed(end,2)];
        acc_mcs = [acc_mcs;b.parameters_trimmed(end,1)];
        str_fts = [str_fts;b.fts_sorted_trimmed];
    end

    % Remove part of the convergence time:
    %temp1 = [acc_data,acc_mcs,tree_size];
    %idx = size(temp1,1) - find(flip(sum(diff(temp1)'))~=0,1) + 1;
    idx = length(time);
    %if isempty(idx)
    %   idx = 1; 
    %end

    r(j).time = time(1:idx);
    r(j).tree_size = tree_size(1:idx);
    r(j).acc_data = acc_data(1:idx);
    r(j).acc_mcs = acc_mcs(1:idx);
    r(j).pop_size = pop_size;
    r(j).str_fts = str_fts(1:idx);
end
end
% Get last results (algorithm output):
function r = get_all_results_output(c)
r = [];
noise = [];
pop_size = [];
acc_before_noise = [];
acc_after_noise = [];
acc_after_cleaning = [];
tree_size = [];
time_convergence = [];
acc_mcs_gt = [];
for j = 1 : length(c)
    a = c(j);
    b = load(strcat(['',a.path,'/output.mat']));
    
    noise = [noise;b.noise];
    pop_size = [pop_size;c(j).pop_size];
    acc_before_noise = [acc_before_noise;b.acc_data_before_noise];
    acc_after_noise = [acc_after_noise;b.acc_data_after_noise];
    %acc_after_cleaning = [acc_after_cleaning;b.acc_data_after_cleaning];
    time_convergence = [time_convergence;b.conv_time];
    acc_mcs_gt = [acc_mcs_gt;b.acc_mcs];
    tree_size = [tree_size;b.tree_size];
    
end

r.noise = noise;
r.pop_size = pop_size;
r.acc_before_noise = acc_before_noise;
r.acc_after_noise = acc_after_noise;
%r.acc_after_cleaning = acc_after_cleaning;
r.time_convergence = time_convergence;
r.acc_mcs_gt = acc_mcs_gt;
r.tree_size = tree_size;


end
function r = get_all_all_results(c)
for j = 1 : length(c)
    a = c(j);
    pop_size  = a.pop_size;
    acc_data  = nan(a.pop_size,length(a.files));
    tree_size = nan(a.pop_size,length(a.files));
    time      = nan(a.pop_size,length(a.files));
    acc_mcs   = nan(a.pop_size,length(a.files));
    str_fts   = {};
    for i = 1 : length(a.files)
        b = load(strcat(['',a.path,'/',a.files{i},'']));
        acc_data(:,i)  = b.parameters_trimmed(:,3);
        time(:,i)      = b.time;
        tree_size(:,i) = b.parameters_trimmed(:,2);
        acc_mcs(:,i)   = b.parameters_trimmed(:,1);
        str_fts        = [str_fts,num2cell(b.fts_sorted_trimmed,2)];
    end

    r(j).time = time;
    r(j).tree_size = tree_size;
    r(j).acc_data = acc_data;
    r(j).acc_mcs = acc_mcs;
    r(j).pop_size = pop_size;
    r(j).str_fts = str_fts;
end
end
%%
function createfigure(XMatrix1, YMatrix1, ZMatrix1)
%CREATEFIGURE(XMatrix1, YMatrix1, ZMatrix1)
%  XMATRIX1:  matrix of x data
%  YMATRIX1:  matrix of y data
%  ZMATRIX1:  matrix of z data

%  Auto-generated by MATLAB on 14-Jun-2021 11:54:26

% Create figure
figure1 = figure('PaperUnits','inches',...
    'PaperSize',[4.66666666666667 4.59722222222222],...
    'Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.178034050570963 0.11 0.726965949429037 0.85257348915056]);
hold(axes1,'on');

% Create multiple lines using matrix input to plot3
plot31 = plot3(XMatrix1,YMatrix1,ZMatrix1,'Parent',axes1,'Marker','o',...
    'LineStyle','none');
set(plot31(1),'DisplayName','Gen. 20','MarkerFaceColor',[1 0 0],...
    'Color',[1 0 0]);
set(plot31(2),'DisplayName','Gen. 10','MarkerFaceColor',[0 1 0],...
    'Color',[0 1 0]);
set(plot31(3),'DisplayName','Gen. 3','MarkerFaceColor',[0 0 1],...
    'Color',[0 0 1]);
set(plot31(4),'DisplayName','Gen. 2','MarkerFaceColor',[1 0 1],...
    'Color',[1 0 1]);
set(plot31(5),'DisplayName','Gen. 1','MarkerFaceColor',[0 0 0],...
    'Color',[0 0 0]);

% Create zlabel
zlabel('\lambda_c');

% Create ylabel
ylabel('\lambda_d');

% Create xlabel
xlabel('\lambda_s');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 100]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.401363326516701 1.0013633265167]);
view(axes1,[-44.9027700831025 50.3333309690156]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',16,'XTick',...
    [0 20 40 60 80 100],'YTick',[0.5 0.6 0.7 0.8 0.9 1],'ZTick',...
    [0 0.2 0.4 0.6 0.8 1]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.569915707893895 0.411144274733941 0.258928571428571 0.274924471299094]);

end