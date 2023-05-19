function visualization_ft_graph(tree_string,colors_BEs)
%clc, clear all, close all

if ~exist('colors_BEs','var')
    colors_BEs = [];
end
%tree_string = 'OR(BE6,BE13,BE4,BE2,BE9,BE17,BE18,BE20,AND(OR(BE3,BE11,BE10),OR(BE12,BE19)),OR(BE13,BE20,BE4,BE2,BE17,BE20,BE4,AND(BE12,OR(BE22,BE11)),BE6,BE9,BE2,BE18,BE13,AND(OR(BE15,BE16,BE14),OR(BE5,BE21,BE23)),BE4,BE17,BE6,AND(OR(BE3,BE19,BE11,BE10),OR(BE12,BE22)),AND(OR(BE5,BE21),OR(BE15,BE16,BE23)),BE18,BE13,BE9,BE2,BE20))';
%colors_BEs([2 4 6 17 9 18 20 13 4]) = 'r';
%colors_BEs([1 3 5 25 19 21 26 8 10 12 14 16 7 11 23 24 22 15]) = 'b';


% Identify the total number of BEs in this FT:
BE_idx = strfind(tree_string,'BE');
c_idx = strfind(tree_string,',');
p_idx = strfind(tree_string,')');
BEs_all = {};
for i = 1 : length(BE_idx)
    c1 = p_idx(find(p_idx>BE_idx(i),1));
    c2 = c_idx(find(c_idx>BE_idx(i),1));
    c = min([c1 c2]);
    BEs_all = [BEs_all,tree_string(BE_idx(i):c-1)];
end

BEs_all = unique(BEs_all)';
%BEs_all = BEs_all';

if ~isempty(colors_BEs)
    % Associate colors to each BE:
    color_per_BE = {};
    for i = 1 : length(BEs_all)
        color_per_BE = [color_per_BE;colors_BEs(str2num(BEs_all{i}(3:end)))];
    end
end
%%  Identify parenthesis:
tree_string_copy = tree_string;
ft_structure = {};
cont_and = length(strfind(tree_string_copy,'AND('))+1;
cont_or = length(strfind(tree_string_copy,'OR('))+1;
list_elements_ft = {};
while ~isempty(i)
    clearvars tree_string_copy_ele
    i = strfind(tree_string_copy,'(');j = strfind(tree_string_copy,')');
    k = j-i';
    k(k<0) = Inf;
    [a,b]= find(k == min(k(:)));
    for l = 1 : length(a)
        p = [i(a(l)),j(b(l))];
    end
    
    for l = 1 : size(p,1)
        all_AND_gates = strfind(tree_string_copy,'AND(');
        all_OR_gates = strfind(tree_string_copy,'OR(');
        % Identify the gate:
        c_and = abs(all_AND_gates-p(l,1));
        c_or = abs(all_OR_gates-p(l,1));
        
        switch_idx = [];
        if isempty(c_and)
            switch_idx = 0;
        elseif isempty(c_or)
            switch_idx = 1;
        elseif min(c_and)<min(c_or) 
            switch_idx = 1;
        elseif min(c_and) > min(c_or) 
            switch_idx = 0;
        end
        
        if switch_idx == 1
            cont_and = cont_and-1;
            % It is an AND gate:
            [~,a]= min(c_and);
            tree_string_copy_ele = extract_elements(tree_string_copy(p(l,1):p(l,2)));
            %ft_structure = [ft_structure; {strcat(['AND_',num2str(cont_and),''])} , {tree_string_copy(p(l,1):p(l,2))}  ];
            ft_structure = [ft_structure; {strcat(['AND_',num2str(cont_and),''])} , {tree_string_copy_ele}];
            tree_string_copy = replaceBetween(tree_string_copy,all_AND_gates(a),p(l,2),strcat(['AND_',num2str(cont_and),'']));
            list_elements_ft = [list_elements_ft;strcat(['AND_',num2str(cont_and),''])];
        elseif switch_idx == 0
            % It is an OR gate:
            cont_or = cont_or-1;
            % It is an AND gate:
            [~,a]= min(c_or);
            tree_string_copy_ele = extract_elements(tree_string_copy(p(l,1):p(l,2)));
            %ft_structure = [ft_structure; {strcat(['OR_',num2str(cont_or),''])} , {tree_string_copy(p(l,1):p(l,2))}  ];
            ft_structure = [ft_structure; {strcat(['OR_',num2str(cont_or),''])} , {tree_string_copy_ele}];
            tree_string_copy = replaceBetween(tree_string_copy,all_OR_gates(a),p(l,2),strcat(['OR_',num2str(cont_or),'']));
            list_elements_ft = [list_elements_ft;strcat(['OR_',num2str(cont_or),''])];
        else
            error('Check veriables c_and and c_or.')
        end
    end
    i = strfind(tree_string_copy,'(');
end
galileo_ft = flip(ft_structure,1);
list_elements_ft = unique(list_elements_ft);

%%
list_ele = [list_elements_ft;BEs_all];
list_ele_coor = zeros(length(list_ele),length(list_ele));
for i = 1 : size(galileo_ft,1)
    idx1 = find(strcmp(list_ele, galileo_ft{i,1}));
    for j = 1 : length(list_ele)
        if any(cell2mat(cellfun(@(x) strcmp(x,list_ele{j}),galileo_ft{i,2},'UniformOutput',false)))
            list_ele_coor(j,idx1) = 1;
        end
    end
end
clearvars temp1 A

a = [];cont = 0;
list_ele_coor2 = ones(size(list_ele_coor)+1)*Inf;
list_ele_coor2(1:size(list_ele_coor,1),1:size(list_ele_coor,1)) = list_ele_coor;
list_ele_coor2(end,1:size(list_ele_coor,1)) = 1:size(list_ele_coor,1);
list_ele_coor2(1:size(list_ele_coor,1),end) = 1:size(list_ele_coor,1);
while ~isempty(list_ele_coor2)
    j = find(sum(list_ele_coor2(1:end-1,1:end-1)') == 0);
    %j = list_ele_coor2(end,j);
    for i = j
        parent = list_ele_coor2(end,i);
        children = list_ele_coor2(find(list_ele_coor2(1:end-1,i)),end);
        if ~isempty(children)
            if isempty(a)
                a = [ones(length(children),1) parent*ones(length(children),1) [2:length(ones(length(children)+1,1))]' children];
            else
                a = [a;a(find(parent==a(:,4)),3)*ones(length(children),1) parent*ones(length(children),1) [1:length(ones(length(children),1))]'+max(a(:,3)) children];
            end
        end
    end
    list_ele_coor2(:,j) = [];list_ele_coor2(j,:) = [];
end

names = {};
b = unique([a(:,1);a(:,3)]);
for i = 1 : length(b)
    idx = unique(a(find(b(i)==a(:,1)),2));
    if ~isempty(idx)
        list_ele{idx}(strfind(list_ele{idx},'_')) = ' ';
        names = [names;list_ele{idx}];
    else
        idx = unique(a(find(b(i)==a(:,3)),4));
        list_ele{idx}(strfind(list_ele{idx},'_')) = ' ';
        names = [names;list_ele{idx}];
    end
    
end

%% Plot the tree shape graph
try
    % Make the second type of graph model:
    type2ftmodel(a,names,colors_BEs,BEs_all)
catch
    % Make the first type of graph model:
    type1ftmodel(a,names,colors_BEs,BEs_all)
end

end
%% Type # 1 FT graph model:
function type1ftmodel(a,names,colors_BEs,BEs_all)
%close all
figure('color',[1 1 1])
s = a(:,1);
t = a(:,3);
G = graph(s,t);
h = plot(G);
nl = names;
h.NodeLabel = '';
xd = get(h, 'XData');
yd = get(h, 'YData');

if ~isempty(colors_BEs)
    for i = 1 : length(nl)
        idx = cell2mat(cellfun(@(x) strcmp(x,nl{i}),BEs_all,'UniformOutput',false));
        if any(idx)
            text(xd(i), yd(i), nl{i}, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle','color',colors_BEs(str2num(BEs_all{idx}(3:end))))
        else
            text(xd(i), yd(i), nl{i}, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
        end
    end
else
    text(xd,yd, nl, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
end
end
%% Type # 2 FT graph model:
function type2ftmodel(a,names,colors_BEs,BEs_all)

x_shift = 0.1;

figure
s = a(:,1);
t = a(:,3);
G = graph(s,t);
h = plot(G);
nl = names;
xd = get(h, 'XData');
yd = get(h, 'YData');
close(gcf)


g = [s t];
p = [xd(2:end)' yd(2:end)'];

% Identify the intermediate highs:
if length(unique(p(:,2)))>1
    diff_high = mode(diff(unique(p(:,2))))/2;
else
    diff_high = 0.5;
end

figure('color',[1 1 1],'position',[1749 540 1244 695]), hold on
% Plot all the elements:
plot(xd,yd,'marker','.','markersize',16,'LineStyle','none','color','k')

% plot the horizontal & vertical lines
for i = [unique(s)]'
    
    idx = find(i==g(:,1));
    t = unique(p(idx,2)) + diff_high;
    plot(p(idx,1),ones(length(idx),1)*t,'k-')
    % Top:
    if i == 1
        plot(ones(1,2)*xd(1),[t (unique(p(idx,2)) + diff_high*2)],'k-')
    else
        plot(ones(1,2)* p(unique(g(idx,1)) == g(:,2),1) ,[t p(unique(g(idx,1)) == g(:,2),2)],'k-')
    end
    
    % Bottom:
    for j = 1:length(idx)
        plot(ones(1,2)*p(idx(j),1),[t p(idx(j),2)],'k-')
    end
    p_xcoor = [g(idx,2) p(idx,1)];
    
end

% Add the labels:
if ~isempty(colors_BEs)
    for i = 1 : length(nl)
        idx = cell2mat(cellfun(@(x) strcmp(x,nl{i}),BEs_all,'UniformOutput',false));
        if any(idx)
            text(xd(i)+x_shift, yd(i), nl{i}, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle','color',colors_BEs(str2num(BEs_all{idx}(3:end))))
            plot(xd(i),yd(i),'color',colors_BEs(str2num(BEs_all{idx}(3:end))) ,'marker','.','markersize',16 )
        else
            text(xd(i)+x_shift, yd(i), nl{i}, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
        end
    end
else
    text(xd+x_shift,yd, nl, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
end
set(gca,'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1])

% Add the top event:
plot(xd(1),max(yd)+diff_high,'marker','.','markersize',16,'color','k')
plot(ones(1,2)*xd(1),[yd(1) yd(1)+diff_high],'marker','.','markersize',16,'color','k')
text(xd(1)+x_shift,max(yd)+diff_high, 'TOP EVENT', 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')

end
%% Extract elements from text
function b = extract_elements(a)
c = sort([strfind(a,'(') strfind(a,')') strfind(a,',')]);
b = {};
for i = 1 : length(c) - 1
    b = [b,a(c(i)+1:c(i+1)-1)];
end
end