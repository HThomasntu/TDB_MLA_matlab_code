% Copyright, M.Bencsik, H.Thomas, 2023

load polygonal_DF_areas.mat
load 27.07.21_sample_2_500hz.mat

extra_bee_X = [0.22595 0.35176 0.43393 0.51662 0.67263 0.84951 0.95725 1.0283 1.1045 1.1858 1.2203 1.2591 1.1261 0.92683]';
extra_bee_Y = [2.2138 2.1697 1.9276 1.8527 1.6698 1.5123 1.5476 1.2822 1.1823 0.99052 0.9708 0.86873 0.78912 0.70083]';


 x = mite_X;
 y = mite_Y;
 
 xx = [bee_X; extra_bee_X];
 yy = [bee_Y; extra_bee_Y];
 
 xxx = bg_X;
 yyy = bg_Y;


% MITE DATA:

mite_data=[x,y];
    plot(mite_data(:,1),mite_data(:,2),'*')
       shrink_f = 0.1;   % shrink factor to tighten the boundary. if increase it will trace outer points exactly and tighten the boundary. 
 
       jb = boundary(mite_data(:,1),mite_data(:,2),shrink_f);
hold on;
        plot(mite_data(jb,1),mite_data(jb,2));
        fill(mite_data(jb,1),mite_data(jb,2),'r');
  hold on;
         plot(df_x,df_y,'o', 'MarkerEdgeColor', 'green');
 hold on 
        plot(x,y,'o','MarkerEdgeColor','black')

%%%-- check any point come inside boundary --%%
[m_bound, n_bound]=size(jb);
%{{
for i=1:1:m_bound
    x_bp(i,1)=mite_data(jb(i,1),1);
    y_bp(i,1)=mite_data(jb(i,1),2);
end
%}
    bound_cord=[x_bp y_bp];
    ans_logic_mite=inpolygon(df_x,df_y,x_bp,y_bp);

if ans_logic_mite==1
    disp(" The point falls inside the polygon");
else
    disp(" The point falls outside the polygon");
end


% BEE DATA: 

bee_data=[xx,yy];
    plot(bee_data(:,1),bee_data(:,2),'*')
       shrink_f = 0.1;   % shrink factor to tighten the boundary. if increase it will trace outer points exactly and tighten the boundary. 
 
       jb2 = boundary(bee_data(:,1),bee_data(:,2),shrink_f);
hold on;
        plot(bee_data(jb2,1),bee_data(jb2,2));
        fill(bee_data(jb2,1),bee_data(jb2,2),'k');
  hold on;
         plot(df_x,df_y,'o', 'MarkerEdgeColor', 'green');
         hold on 
        plot(xx,yy,'o','MarkerEdgeColor','white')
        plot(xx,yy,'.')

%%%-- check any point come inside boundary --%%
[m_bound2, n_bound2]=size(jb2);
%{{
for i2=1:1:m_bound2
    xx_bp(i2,1)=bee_data(jb2(i2,1),1);
    yy_bp(i2,1)=bee_data(jb2(i2,1),2);
end
%}
    bound_cord2=[xx_bp yy_bp];
    ans_logic_bee=inpolygon(df_x,df_y,xx_bp,yy_bp);

if ans_logic_bee==1
    disp(" The point falls inside the polygon");
else
    disp(" The point falls outside the polygon");
end


% BG DATA: 

bg_data=[xxx,yyy];
    plot(bg_data(:,1),bg_data(:,2),'*')
       shrink_f = 0.1;   % shrink factor to tighten the boundary. if increase it will trace outer points exactly and tighten the boundary. 
 
       jb3 = boundary(bg_data(:,1),bg_data(:,2),shrink_f);
hold on;
        plot(bg_data(jb3,1),bg_data(jb3,2));
        fill(bg_data(jb3,1),bg_data(jb3,2),'b');
  hold on;
         plot(df_x,df_y,'o', 'MarkerEdgeColor', 'green');
         hold on 
        plot(xxx,yyy,'o','MarkerEdgeColor','yellow')

%%%-- check any point come inside boundary --%%
[m_bound3, n_bound3]=size(jb3);
%{{
for i3=1:1:m_bound3
    xxx_bp(i3,1)=bg_data(jb3(i3,1),1);
    yyy_bp(i3,1)=bg_data(jb3(i3,1),2);
end
%}
    bound_cord3=[xxx_bp yyy_bp];
    ans_logic_bg=inpolygon(df_x,df_y,xxx_bp,yyy_bp);

if ans_logic_bg==1
    disp(" The point falls inside the polygon");
else
    disp(" The point falls outside the polygon");
    
end

% which masked area does each point fall into?

new_matrix = zeros(3, length(df_x));

for i = 1:length(df_x)
    ans_logic_mite = inpolygon(df_x(1,i),df_y(1,i),x_bp,y_bp);
    ans_logic_bee = inpolygon(df_x(1,i),df_y(1,i),xx_bp,yy_bp);
    ans_logic_bg = inpolygon(df_x(1,i),df_y(1,i),xxx_bp,yyy_bp);
    if ans_logic_mite == 1 && ans_logic_bee == 0 && ans_logic_bg == 0
        disp ("the point falls inside mite")
        rule_all(1,i) = 999;
        new_matrix(1,i)= df_x(1,i);
        new_matrix(2,i)= df_y(1,i);
        new_matrix(3,i)= rule_all(1,i);
    elseif ans_logic_mite == 0 && ans_logic_bee == 1 && ans_logic_bg == 0
        disp ("the point falls inside bee")
        rule_all(1,i) = 998;
        new_matrix(1,i)= df_x(1,i);
        new_matrix(2,i)= df_y(1,i);
        new_matrix(3,i)= rule_all(1,i);
    elseif ans_logic_mite == 0 && ans_logic_bee == 0 && ans_logic_bg == 1
        disp ("the point falls inside background ")
        rule_all(1,i) = 997;
        new_matrix(1,i)= df_x(1,i);
        new_matrix(2,i)= df_y(1,i);
        new_matrix(3,i)= rule_all(1,i);
    else
        disp ("the point falls in no category")
        rule_all(1,i) = 996;
        new_matrix(1,i)= df_x(1,i);
        new_matrix(2,i)= df_y(1,i);
        new_matrix(3,i)= rule_all(1,i); 
    end
end

new_matrix;

% view points as 2dfts to check whether they have fallen in the correct
% area:

[a b] = find(new_matrix(3,:) == 996); % no category
[c d] = find(new_matrix(3,:) == 997); % background
[e f] = find(new_matrix(3,:) == 998); % bee
[g h] = find(new_matrix(3,:) == 999); % mite
