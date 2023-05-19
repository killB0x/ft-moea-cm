clc, clear all, close all

%% Initial parameters:
%cases_name = {'container_seal_design_noise_0per','covid-19_noise_0per','monopropellant_propulsion_system_noise_0per','presure_tank_noise_0per'};
cases_name = {'monopropellant_propulsion_system_noise_0per'};
gt_all = 23;%[10,33,23,11];
save_plot = 0;

for n = 1 : length(cases_name)
    save_plots_dir = strcat(['',cases_name{n},'']);
    if ~isfolder(save_plots_dir)
        mkdir(save_plots_dir)
    end
    c = {'_MCS1TS1Acc1','_MCS0TS0Acc1'};
    check_folders = {};
    for i = 1:length(c)
        if i == 1
            a = 'sdc';
        else
            a = 's';
        end
        check_folders = [check_folders;strcat(['data_for_figures/fig6/',a,'/',cases_name{n},'',c{i},'_sfv0'])];
    end
    
    compare = {};
    output_results = {};
    for i = 1 : length(check_folders)
        output_results = [output_results;get_all_results_output(get_dir_files_str(check_folders{i}))];
        compare = [compare;get_all_results(get_dir_files_str(check_folders{i}))];
    end

    m = {'sdc','d'};
    gt = gt_all(n);% Ground truth fault tree size 
    
    
    plot_setup.fig_pos        = [1477         241         213         197];
    plot_setup.axis_pos       = [0.3141    0.3421    0.6144    0.6325];
    plot_setup.text_size      = 16;
    plot_setup.legend_pos     = [0.2142    0.7810    0.7230    0.2059];
    plot_setup.markers        = 'oxs^+';
    plot_setup.font_name      = 'Arial Unicode MS';
    plot_setup.save_plots_dir = save_plots_dir;
    
    %%
    % Accuracy Cut Sets vs Population Size (with noise):
    acc_cut_sets_vs_pop_size_with_noise_hist(output_results,m,save_plot,plot_setup)
    % Accuracy Dataset vs Population Size (with noise):
    acc_dataset_vs_pop_size_with_noise(output_results,m,save_plot,plot_setup)
    % Tree Size vs Accuracy Cut Sets (with noise):
    tree_size_vs_acc_mcs_with_noise(output_results,m,save_plot,plot_setup,gt)
    % Tree Size vs Accuracy dataset (with noise):
    tree_size_vs_acc_dataset_with_noise(output_results,m,save_plot,plot_setup,gt)
    % Tree Size vs Population Size:
    tree_size_vs_pop_size_with_noise(output_results,m,save_plot,plot_setup,gt)
    % Time vs Population Size:
    time_vs_pop_size(compare,m,save_plot,plot_setup)
    
end
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
    temp1 = [acc_data,acc_mcs,tree_size];
    idx = size(temp1,1) - find(flip(sum(diff(temp1)'))~=0,1) + 1;
    
    if isempty(idx)
       idx = 1; 
    end

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
    acc_after_cleaning = [acc_after_cleaning;b.acc_data_after_cleaning];
    time_convergence = [time_convergence;b.conv_time];
    acc_mcs_gt = [acc_mcs_gt;b.acc_mcs];
    tree_size = [tree_size;b.tree_size];
    
end

r.noise = noise;
r.pop_size = pop_size;
r.acc_before_noise = acc_before_noise;
r.acc_after_noise = acc_after_noise;
r.acc_after_cleaning = acc_after_cleaning;
r.time_convergence = time_convergence;
r.acc_mcs_gt = acc_mcs_gt;
r.tree_size = tree_size;


end
% Accuracy Cut Sets vs Population Size:
function acc_cut_sets_vs_pop_size_with_noise_hist(compare,m,save_plot,plot_setup)
fig_pos = plot_setup.fig_pos;
axis_pos = plot_setup.axis_pos;
text_size = plot_setup.text_size;
font_name=plot_setup.font_name;

a = [];
for i = 1 : length(compare)
    a = [a,[compare{i}.pop_size]];
end

b=[];d = {};
for i = 1 : length(compare)
    b = [b,compare{i}.acc_mcs_gt];
    d = [d,compare{i}.acc_mcs_gt];
end

e = {};
for i = 1 : length(d)
    e = [e;num2cell(repmat(m{i},length(d{i}),1),2)];
end

a = reshape(a,[],1);b = reshape(b,[],1);
x = categorical(a,unique(a));
figure('color',[1 1 1],'position',fig_pos),hold on
boxchart(x,1-b,'GroupByColor',e)
xlabel('Population size')
ylabel('\phi_c')
%ylim([min(b) 1])
box on
set(gca,'FontName',font_name,'FontSize',text_size,'position',axis_pos)
leg1 = legend('show');
set(leg1,'location','best')
%set(leg1,'position',plot_setup.legend_pos)
ylim([-0.01 max(1-b)])

if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',plot_setup.save_plots_dir ,'/acc_cut_sets_vs_pop_size']),'-dpdf','-r0')
    close(gcf)
end
end
% Accuracy Dataset vs Population Size (with noise):
function acc_dataset_vs_pop_size_with_noise(compare,m,save_plot,plot_setup)
fig_pos = plot_setup.fig_pos;
axis_pos = plot_setup.axis_pos;
text_size = plot_setup.text_size;
font_name=plot_setup.font_name;

a = [];
for i = 1 : length(compare)
    a = [a,[compare{i}.pop_size]];
end

b=[];d = {};
for i = 1 : length(compare)
    b = [b,compare{i}.acc_after_noise];
    d = [d,compare{i}.acc_after_noise];
end

e = {};
for i = 1 : length(d)
    e = [e;num2cell(repmat(m{i},length(d{i}),1),2)];
end

a = reshape(a,[],1);b = reshape(b,[],1);
x = categorical(a,unique(a));
figure('color',[1 1 1],'position',fig_pos),hold on
boxchart(x,1-b,'GroupByColor',e)
xlabel('Population size')
ylabel('\phi_d')
%ylim([min(b) 1])
box on
set(gca,'FontName',font_name,'FontSize',text_size,'position',axis_pos)
leg1 = legend('show');
set(leg1,'location','best')
ylim([-0.01 max(1-b)])
%set(leg1,'position',plot_setup.legend_pos)

if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',plot_setup.save_plots_dir ,'/acc_data_vs_pop_size']),'-dpdf','-r0')
    close(gcf)
end
end
% Tree Size vs Accuracy Cut Sets (with noise):
function tree_size_vs_acc_mcs_with_noise(compare,m,save_plot,plot_setup,gt)
fig_pos = plot_setup.fig_pos;
axis_pos = plot_setup.axis_pos;
text_size = plot_setup.text_size;
font_name=plot_setup.font_name;

a = [];
for i = 1 : length(compare)
    a = [a,[compare{i}.pop_size]];
end

figure('color',[1 1 1],'position',fig_pos),hold on
b=[];
compare = flip(compare);m = flip(m);
for i = 1 : length(compare)
    c = [double(compare{i}.tree_size),compare{i}.acc_mcs_gt];
    plot(1-c(:,2),c(:,1),'MarkerSize',8,'DisplayName',m{i},'Marker',plot_setup.markers(i),'LineStyle','none')
    b = [b;c];
end

for i = 1 : length(gt)
    % Plot the ground truth:
    %plot([min(b(:,2)) 1],gt(i)*ones(2,1),'r-.')
    plot([0 max(get(gca,'XTick'))],gt(i)*ones(2,1),'r-.')
end

xlabel('\phi_c')
ylabel('\phi_s')
%xlim([min(b(:,2)) 1])
%xlim([0 1])
box on
set(gca,'FontName',font_name,'FontSize',text_size,'position',axis_pos,'YTick',[0:20:60])
leg1 = legend(m);
set(leg1,'location','best')
%set(leg1,'position',plot_setup.legend_pos)


if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',plot_setup.save_plots_dir ,'/tree_size_vs_acc_mcs']),'-dpdf','-r0')
    close(gcf)
end
end
% Tree Size vs Accuracy dataset (with noise):
function tree_size_vs_acc_dataset_with_noise(compare,m,save_plot,plot_setup,gt)
fig_pos = plot_setup.fig_pos;
axis_pos = plot_setup.axis_pos;
text_size = plot_setup.text_size;
font_name=plot_setup.font_name;

a = [];
for i = 1 : length(compare)
    a = [a,[compare{i}.pop_size]];
end

figure('color',[1 1 1],'position',fig_pos),hold on
b=[];

compare = flip(compare);m = flip(m);
for i = 1 : length(compare)
    c = [double(compare{i}.tree_size) compare{i}.acc_after_noise];
    plot(1-c(:,2),c(:,1),'MarkerSize',8,'DisplayName',m{i},'Marker',plot_setup.markers(i),'LineStyle','none')
    b = [b;c];
end

for i = 1 : length(gt)
    % Plot the ground truth:
    %plot([min(b(:,2)) 1],gt(i)*ones(2,1),'r-.')
    plot([0 max(get(gca,'XTick'))],gt(i)*ones(2,1),'r-.')
end

xlabel('\phi_d')
ylabel('\phi_s')
%xlim([min(b(:,2)) 1])
%xlim([0 1])
box on
set(gca,'FontName',font_name,'FontSize',text_size,'position',axis_pos,'YTick',[0:20:60])
leg1 = legend(m);
set(leg1,'location','best')
%set(leg1,'position',plot_setup.legend_pos)


if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',plot_setup.save_plots_dir ,'/tree_size_vs_acc_dataset']),'-dpdf','-r0')
    close(gcf)
end
end
% Tree Size vs Population Size (with noise):
function tree_size_vs_pop_size_with_noise(compare,m,save_plot,plot_setup,gt)
fig_pos = plot_setup.fig_pos;
axis_pos = plot_setup.axis_pos;
text_size = plot_setup.text_size;
font_name=plot_setup.font_name;

a = [];
for i = 1 : length(compare)
    a = [a,[compare{i}.pop_size]];
end

b=[];d = {};
for i = 1 : length(compare)
    b = [b,double(compare{i}.tree_size)];
    d = [d,double(compare{i}.tree_size)];
end

e = {};
for i = 1 : length(d)
    e = [e;num2cell(repmat(m{i},length(d{i}),1),2)];
end

a = reshape(a,[],1);b = reshape(b,[],1);
x = categorical(a,unique(a));
figure('color',[1 1 1],'position',fig_pos),hold on
boxchart(x,b,'GroupByColor',e)
xlabel('Population size')
ylabel('\phi_s')
box on
set(gca,'FontName',font_name,'FontSize',text_size,'position',axis_pos,'YTick',[0:20:60])


for i = 1 : length(gt)
    % Plot the ground truth:
    plot([0 max(a)],gt(i)*ones(2,1),'r-')
end


%leg1 = legend(m);
leg1 = legend('show');
set(leg1,'location','best','String',leg1.String(1:end-1))
%set(leg1,'position',plot_setup.legend_pos)

if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',plot_setup.save_plots_dir ,'/tree_size_vs_pop_size']),'-dpdf','-r0')
    close(gcf)
end

end
% Time vs Population Size:
function time_vs_pop_size(compare,m,save_plot,plot_setup)
fig_pos = plot_setup.fig_pos;
axis_pos = plot_setup.axis_pos;
text_size = plot_setup.text_size;
font_name=plot_setup.font_name;

a = [];
for i = 1 : length(compare)
    a = [a,[compare{i}.pop_size]];
end

b=[];d = {};
for i = 1 : length(compare)
    c=[];
    for j = 1 : length(compare{i})
        
        % Check the time that took for convergence:
        c = [c,compare{i}(j).time(end)/60];
    end
    b = [b,c];
    d = [d,c];
end

e = {};
for i = 1 : length(d)
    e = [e;num2cell(repmat(m{i},length(d{i}),1),2)];
end

x = categorical(a,unique(a));
figure('color',[1 1 1],'position',fig_pos),hold on
boxchart(x,b,'GroupByColor',e)
xlabel('Population size')
ylabel('Time (min)')
box on
set(gca,'FontName',font_name,'FontSize',text_size,'position',axis_pos)
leg1 = legend('show');
set(leg1,'location','best')
%set(leg1,'position',plot_setup.legend_pos)
if save_plot == 1
    h = gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,strcat(['',plot_setup.save_plots_dir ,'/time_vs_pop_size']),'-dpdf','-r0')
    close(gcf)
end
end