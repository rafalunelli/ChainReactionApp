clear all;
close all;
clc;


gridSize = [1000,1000];

%quantidades iniciais
N235 = 100;
N238 = 50;
Nn = 1;
%raios: U235,U238,nêutron,Bário,Criptônio
radius = [15,15,5,10,10]; 

%nêutron
v0_n = [300,0]; 

%Criptônio-91
v0_Kr = [150,150]; 
%Bário-142
v0_Ba = [150,150]; 

objects = [];
tem_nan = false;

%create uranium
while size(objects,1)<(N238+N235)
    tem_nan = false;
    
    if size(objects,1) < N235
        new_obj = {'U235',rand*gridSize(1),rand*gridSize(2),0,0};
        current = 1;
    else
        new_obj = {'U238',rand*gridSize(1),rand*gridSize(2),0,0};
        current = 2;
    end
    
    if ~isempty(objects)
        for obj_i = 1:size(objects,1)
            if isequal(objects{obj_i,1},'U235')
                last = 1; 
            else
                last = 2;
            end
            [xout,yout] = circcirc(new_obj{2},new_obj{3},radius(current),objects{obj_i,2},objects{obj_i,3},radius(last));
            if ~isnan(xout)
                tem_nan = true;
                break;                
            end
        end
    end    
    if (~tem_nan)
        objects = [objects;new_obj];
    end
end

%create neutrons
for i = 1:Nn
   s0 = [-300,gridSize(2)/2];
   v0x = v0_n.*rand(1,2);
   new_obj = {'n',s0(1),s0(2),v0x(1),v0x(2)};
   
   objects(end+1,:)=new_obj; 
end

ax = gca;
hold on
plot_objs = [];
for obj_id = 1:size(objects)
        if isequal(objects{obj_id,1},'U235')
            p= viscircles([objects{obj_id,2},objects{obj_id,3}],radius(1),'Color','r');
            %text(objects{obj_id,2},objects{obj_id,3},num2str(obj_id));
        elseif isequal(objects{obj_id,1},'U238')
            p= viscircles([objects{obj_id,2},objects{obj_id,3}],radius(2),'Color','k');
            %text(objects{obj_id,2},objects{obj_id,3},num2str(obj_id));
        elseif isequal(objects{obj_id,1},'n')
            p = plot(objects{obj_id,2},objects{obj_id,3},'.k','MarkerSize',10);
        end
        plot_objs = [plot_objs;{objects{obj_id,1},p}];
end
axis equal
axis manual;

t0 = 0;
tf = 50;
dt = 0.01;

t = 0;

[TimeNextCollision,NNextCollision,UnextCollision] = nextEvent(objects,radius);
NextEventTime = t + TimeNextCollision;  
while t<tf
    tic
    t=t+dt;
    
    [objects,plot_objs] = UpdatePositions(objects,plot_objs,dt);
    
    if (t>=NextEventTime-dt)
            U235Position = [objects{UnextCollision,2:3}];
            delete(plot_objs{NNextCollision,2});
            delete(plot_objs{UnextCollision,2});            
            objects([NNextCollision,UnextCollision],:)=[];
            plot_objs([NNextCollision,UnextCollision],:)=[];
            
            %criar novos nêutrons:        
            [objects,plot_objs] = create3Neutrons(ax,v0_n,U235Position,objects,plot_objs);
            %criar Ba, Kr
            [objects,plot_objs] = createKrBa(ax,v0_Kr,v0_Ba,U235Position,radius,objects,plot_objs);   
            
            [objects,plot_objs] = RemoveElementsOutside(ax,radius,objects,plot_objs);
            
            [TimeNextCollision,NNextCollision,UnextCollision] = nextEvent(objects,radius);
            NextEventTime = t + TimeNextCollision;  
    end    
    run_time = toc;
    pause(max(1e-10,dt-run_time));
end

function [objects,plot_objs] = RemoveElementsOutside(ax,radius,objects,plot_objs)

objects_positions = [[objects{:,2}]',[objects{:,3}]'];

%raios: U235,U238,nêutron,Bário,Criptônio

XLIM = ax.XLim;
YLIM = ax.YLim;

obj_names = {'n','Kr91','Ba142'};
obj_raios = radius(3:5);

remove_total = [];
for obj_id = 1:3
    
    obj_indices = find(ismember(plot_objs(:,1),obj_names{obj_id}));
    remove_ind = or((objects_positions(obj_indices,1)+obj_raios(obj_id))<XLIM(1),...
        or((objects_positions(obj_indices,1)-obj_raios(obj_id))>XLIM(2),...
        or((objects_positions(obj_indices,2)+obj_raios(obj_id))<YLIM(1),...
        (objects_positions(obj_indices,2)-obj_raios(obj_id))>YLIM(2))));
    
    remove_total = [remove_total;obj_indices(logical(remove_ind))];
end
if ~isempty(remove_total)
    objects(remove_total,:) = [];
    delete([plot_objs{remove_total,2}]);
    plot_objs(remove_total,:) = [];
end


end

function [objects,plot_objs] = UpdatePositions(objects,plot_objs,dt)

    uranios = find(ismember(plot_objs(:,1),'n'));
    subnucleos = find(ismember(plot_objs(:,1),{'Kr91','Ba142'}));
    
    for i = 1:length(uranios)
        pi = uranios(i);
        plot_objs{pi,2}.XData = plot_objs{pi,2}.XData + objects{pi,4}*dt;
        plot_objs{pi,2}.YData = plot_objs{pi,2}.YData + objects{pi,5}*dt;
        objects{pi,2} = plot_objs{pi,2}.XData;
        objects{pi,3} = plot_objs{pi,2}.YData;
    end
    for i = 1:length(subnucleos)
        pi = subnucleos(i);
        
        plot_objs{pi,2}.Children(1).XData = plot_objs{pi,2}.Children(1).XData + objects{pi,4}*dt;
        plot_objs{pi,2}.Children(1).YData = plot_objs{pi,2}.Children(1).YData + objects{pi,5}*dt;
        plot_objs{pi,2}.Children(2).XData = plot_objs{pi,2}.Children(1).XData;
        plot_objs{pi,2}.Children(2).YData = plot_objs{pi,2}.Children(1).YData;
        objects{pi,2} = objects{pi,2} + objects{pi,4}*dt;
        objects{pi,3} = objects{pi,3} + objects{pi,5}*dt;
    end   
    
end

function [objects,plot_objs] = createKrBa(ax,v0_Kr,v0_Ba,s0,radius,objects,plot_objs)
axes(ax);
%Bário
v0xBa = v0_Ba(1)*(2*randi(2)-3);
v0yBa = v0_Ba(2)*(2*randi(2)-3);
new_Ba = {'Ba142',s0(1),s0(2),v0xBa,v0yBa};
p1 = viscircles([s0(1),s0(2)],radius(4),'Color','b');

objects(end+1,:) = new_Ba;
plot_objs = [plot_objs;{'Ba142',p1}];


%Kriptônio
v0xKr = v0_Kr(1)*(2*randi(2)-3);
v0yKr = v0_Kr(2)*(2*randi(2)-3);
new_Kr = {'Kr91',s0(1),s0(2),v0xKr,v0yKr};
p2 = viscircles([s0(1),s0(2)],radius(5),'Color','g');

objects(end+1,:) = new_Kr;
plot_objs = [plot_objs;{'Kr91',p2}];



end

function [objects,plot_objs] = create3Neutrons(ax,v0,s0,objects,plot_objs)
    axes(ax);
    for i = 1:3
        % v0x = 2*mean(v0)*(2*rand-1);
        % v0y = 2*mean(v0)*(2*rand-1);
        v0x = (2*mean(v0))*(2*randi(2)-3);
        v0y = (2*mean(v0))*(2*randi(2)-3);
        
        objects(end+1,:) = {'n',s0(1),s0(2),v0x,v0y};
        plot_objs = [plot_objs;{'n',plot(s0(1),s0(2),'.k','MarkerSize',10)}];
    end
end

function [TimeNextCollision,NNextCollision,UnextCollision]= nextEvent(Allobjects,radius)
    NNextCollision = -1;
    TimeNextCollision = inf;
    UnextCollision = -1;

    index_neutrons = find(ismember(Allobjects(:,1),'n'));
    index_u235 = find(ismember(Allobjects(:,1),'U235'));
    for nid = 1:length(index_neutrons)        
        for uid = 1:length(index_u235)
           dx = Allobjects{index_neutrons(nid),2}-Allobjects{index_u235(uid),2};
           dy = Allobjects{index_neutrons(nid),3}-Allobjects{index_u235(uid),3};
            
           dvx = Allobjects{index_neutrons(nid),4};
           dvy = Allobjects{index_neutrons(nid),5};
           dv = dvx*dvx + dvy*dvy;
           
           scal = dvx*dx+dvy*dy;
           ups = (scal*scal)- dv*(dx*dx+dy*dy-4*radius(1)*radius(3));
           if(ups>0)&&(scal<0)
               deltaT = (-scal-sqrt(ups))/(dvx^2+dvy^2);
               if ((deltaT<TimeNextCollision)&&deltaT>0)
                   
                   NNextCollision = index_neutrons(nid);
                   TimeNextCollision = deltaT;
                   UnextCollision = index_u235(uid);
               end
           end
        end
    end
    
    
    
end
