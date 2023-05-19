clc, clear all, close all
%% Acc_MCS, Tree Size, Acc_data over the generations:
plot_setup.save_plots_dir = 'figures_in_paper';
case_name_all = {'MCS1TS1Acc1','MCS0TS0Acc1'};
folder_results_all{1} = 'data_for_figures/fig4_5/MCS1TS1Acc1/c1_run_1_pop_size_400/run_0';
folder_results_all{2} = 'data_for_figures/fig4_5/MCS0TS0Acc1/c1_run_1_pop_size_400/run_0';

figure_pos = [1233         450         213         190];
axis_pos   = [0.3141    0.3216    0.6144    0.6455];%[0.3474 0.3076 0.5983 0.6556];
save_plot  = 1;
pop_size   = 400;

font_name = 'Arial Unicode MS';
text_size= 16;
YTick_all = {[0:5:25],[0:50:150]};

for j = 1 : length(folder_results_all)
    case_name = case_name_all{j};
    folder_results = folder_results_all{j};
    files = dir(folder_results);
    files = {files.name}';
    a = [];
    for i = 1:length(files)
        if files{i}(1) ~= 'g'
            a = [a;i];
        end
    end
    files(a) = [];
    g = [];
    for i = 1 : length(files)
        a = strfind(files{i},'_');
        g = [g;str2num(files{i}(a(end)+1:end-4))];
    end
    [~,x] = sort(g,'ascend');
    s = g(x);
    files = files(x);
    
    acc_data_all = {}; tree_size_all = {}; time_ = {}; fts_str={}; acc_mcs_all= {};
    acc_data = zeros(length(files),pop_size);tree_size = zeros(length(files),pop_size);acc_mcs = zeros(length(files),pop_size);
    for i = 1 : length(files)
        a = load(strcat(['',folder_results,'/',files{i},'']));
        acc_mcs(i,:) = a.parameters_trimmed(:,1)';
        tree_size(i,:) = a.parameters_trimmed(:,2)';
        acc_data(i,:) = a.parameters_trimmed(:,3)';
    end
    
    acc_data_all = [acc_data_all,acc_data];
    tree_size_all = [tree_size_all;tree_size];
    acc_mcs_all = [acc_mcs_all;acc_mcs];
    
    % Tree size
    i = 1;
    figure('color',[1 1 1],'position',figure_pos)
    hold on
    plot_distribution_prctile(0:size(tree_size_all{i},1)-1,tree_size_all{i}','Prctile',[25:25:100],'color',[0     0     0]);
    hold on
    yl = get(gca,'YLim');
    plot([0 40],[23 23],'b-','linewidth',1)
    %plot([20 20],yl,'r-','linewidth',1)
    plot(0:size(acc_mcs_all{i},1)-1,max(tree_size_all{i}'),'Marker','o','MarkerFaceColor','r','MarkerSize',3,'LineStyle','none','Color',[1 0 0])
    hold off
    ylabel('\phi_s')
    set(gca,'FontSize',text_size,'FontName',font_name,'position',axis_pos,'XTick',[0:10:40],'YTick',YTick_all{j})
    %grid on, 
    box on
    ylim(yl)
    %title(strcat(['Tree size for Max. Acc. Case D2']),'FontWeight','normal')
    xlabel('Generation')
    xlim([0 40])
    if save_plot == 1
        exportgraphics(gcf,strcat(['',plot_setup.save_plots_dir ,'',case_name,'_convergence_tree_size_over_generations.png']),'Resolution',300)
        savefig(gcf,strcat(['',plot_setup.save_plots_dir ,'',case_name,'_convergence_tree_size_over_generations.fig']))
        close(gcf)
    end
    
    % Accuracy based on the data
    figure('color',[1 1 1],'position',figure_pos)
    hold on
    plot_distribution_prctile(0:size(acc_data_all{i},1)-1,1-acc_data_all{i}','Prctile',[5:10:100],'color',[0     0     0]);
    hold on
    plot(0:size(acc_mcs_all{i},1)-1,1-max(acc_data_all{1}'),'Marker','o','MarkerFaceColor','r','MarkerSize',3,'LineStyle','none','Color',[1 0 0])
    %plot([20 20],[0 1],'r-','linewidth',1)
    hold off
    ylabel('\phi_d')
    set(gca,'FontSize',text_size,'FontName',font_name,'position',axis_pos,'XTick',[0:10:40],'YTick',[0:0.25:1])
    %grid on, 
    box on
    %title(strcat(['Max. Acc. Case D2']),'FontWeight','normal')
    xlabel('Generation')
    xlim([0 40])
    %ylim([0.6 1])
    ylim([0 1])
    if save_plot == 1
        exportgraphics(gcf,strcat(['',plot_setup.save_plots_dir ,'',case_name,'_convergence_acc_data_over_generations.png']),'Resolution',300)
        savefig(gcf,strcat(['',plot_setup.save_plots_dir ,'',case_name,'_convergence_acc_data_over_generations.fig']))
        close(gcf)
    end
    
    
    % Accuracy based on the MCSs
    figure('color',[1 1 1],'position',figure_pos)
    hold on
    plot_distribution_prctile(0:size(acc_mcs_all{i},1)-1,1-acc_mcs_all{i}','Prctile',[5:10:100],'color',[0     0     0]);
    hold on
    plot(0:size(acc_mcs_all{i},1)-1,1-max(acc_mcs_all{1}'),'Marker','o','MarkerFaceColor','r','MarkerSize',3,'LineStyle','none','Color',[1 0 0])
    %plot([20 20],[0 1],'r-','linewidth',1)
    hold off
    ylabel('\phi_c')
    set(gca,'FontSize',text_size,'FontName',font_name,'position',axis_pos,'XTick',[0:10:40],'YTick',[0:0.25:1])
    %grid on, 
    box on
    %title(strcat(['Max. Acc. Case D2']),'FontWeight','normal')
    xlabel('Generation')
    xlim([0 40])
    ylim([0 1])
    if save_plot == 1
        exportgraphics(gcf,strcat(['',plot_setup.save_plots_dir ,'',case_name,'_convergence_acc_mcs_over_generations.png']),'Resolution',300)
        savefig(gcf,strcat(['',plot_setup.save_plots_dir ,'',case_name,'_convergence_acc_mcs_over_generations.fig']))
        close(gcf)
    end
    
end
%