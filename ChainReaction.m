classdef ChainReaction
        
    properties
        running;
        Ax;
        gridSize;
        
        %quantidades:
        N235;
        N238;
        Nn;
        %raios:
        r235;
        r238;
        r239;
        rn;
        rBa;
        rKr;
        
        %velocidades iniciais:
        v0_n;
        v0_Kr;
        v0_Ba;
        %simulação:
        t0;
        dt;
        tf;
        %objetos:
        objects;
        compute_obj;
        
        %image files
        u235_file;
        u235_alpha;
        u238_file;
        u238_alpha;
        u239_file;
        u239_alpha;
        Kr91_file;
        Kr91_alpha;
        Ba142_file;
        Ba142_alpha;
        
    end
    
    methods
        
        function obj = ChainReaction(Ax)
            obj.running = false;
            obj.Ax = Ax;
            
            [obj.u239_file, ~, obj.u239_alpha] =  imread('images\u239.png');
            [obj.u238_file, ~, obj.u238_alpha] =  imread('images\u238.png');
            [obj.u235_file, ~, obj.u235_alpha] =  imread('images\u235.png');
            [obj.Kr91_file, ~, obj.Kr91_alpha] =  imread('images\Kr91.png');
            [obj.Ba142_file, ~, obj.Ba142_alpha] =  imread('images\Ba142.png');
        end
        
        function [obj,n_real,u_real,m_real,t_real] = ChainSim(obj,t0,tf,dt,app)
            obj.running = true;
            
            n_real = zeros((tf-t0)/dt,1);
            u_real = zeros((tf-t0)/dt,1);
            m_real = zeros((tf-t0)/dt,1); 
            t_real = zeros((tf-t0)/dt,1);
            
            n_real(1) = obj.Nn;
            u_real(1) = obj.N235;
            m_real(1) = obj.N238;
            t_real(1) = t0;
            
            obj.dt = dt;
            obj.tf = tf;
            obj.t0 = t0;
            
            t = t0;
            
            [obj,TimeNextCollision,NNextCollision,UnextCollision] = obj.nextEvent();
            NextEventTime = t + TimeNextCollision;  
            
            %while t<tf
            for t_i = 1:(tf-t0)/dt
                
                if ~app.running
                    break
                end
                
                tic
                t = t+dt;

                obj = obj.UpdatePositions();

                if (t>=(NextEventTime-dt))
                    
                        UName = obj.objects{UnextCollision};
                        UPosition = [obj.objects{UnextCollision,2:3}];
                        deleteded_n = obj.objects(NNextCollision,:);
                        
                        delete(obj.objects{NNextCollision,6});
                        delete(obj.objects{UnextCollision,6});            
                        obj.objects([NNextCollision,UnextCollision],:)=[];
                        obj.compute_obj([NNextCollision,UnextCollision],:)=[];
                        
                        if (isequal(UName,'U235'))
                            %criar novos nêutrons:        
                            obj = obj.create3Neutrons(UPosition);
                            %criar Ba, Kr:
                            obj = obj.createKrBa(UPosition); 
                            
                            n_real(t_i) = n_real(t_i-1) + 2;
                            u_real(t_i) = u_real(t_i-1) - 1;
                            m_real(t_i) = m_real(t_i-1);
                            t_real(t_i) = t; 
                        else
                            obj = obj.createU239(UPosition);
                            n_real(t_i) = n_real(t_i-1) - 1;
                            u_real(t_i) = u_real(t_i-1);
                            m_real(t_i) = m_real(t_i-1) - 1;
                            t_real(t_i) = t; 
                        end
                        
                        obj = obj.RemoveElementsOutside();
                        [obj,TimeNextCollision,NNextCollision,UnextCollision] = obj.nextEvent();
                        NextEventTime = t + TimeNextCollision;  
                        
                else
                    n_real(t_i) = n_real(max(1,t_i-1));
                    u_real(t_i) = u_real(max(1,t_i-1));
                    m_real(t_i) = m_real(max(1,t_i-1));
                    t_real(t_i) = t;
                end            
                app.Gauge.Value = t;
                run_time = toc;
                pause(max(1e-10,dt-run_time));
                
            end
            
            position_first = 1;%find(t_real > (300/obj.v0_n(1)),1);
            t_real = t_real(position_first:t_i-1);
            n_real = n_real(position_first:t_i-1);
            u_real = u_real(position_first:t_i-1);
            m_real = m_real(position_first:t_i-1);
            
            obj.running = false;
        end
        
        function obj = GetParameters(obj,gridSize, N235,N238,Nn, r235,r238,r239,rn,rBa,rKr, v0_n,v0_Kr,v0_Ba)
            
            obj.running = false;
            obj.gridSize = gridSize;
            
            obj.N235 = N235;
            obj.N238 = N238;
            obj.Nn = Nn;
            
            obj.r235 = r235;
            obj.r238 = r238;
            obj.r239 = r239;
            obj.rn = rn;
            obj.rBa = rBa;
            obj.rKr = rKr;
            
            obj.v0_n = v0_n;
            obj.v0_Kr = v0_Kr;
            obj.v0_Ba = v0_Ba;
                       
        end
        
        function obj = GenerateElements(obj)
            obj.objects = [];
            obj.compute_obj = [];
            
            %create uranium 235
            while size(obj.objects,1) < obj.N235    
                tem_nan = false;
                s0 = rand(1,2).*obj.gridSize; 

                if ~isempty(obj.objects)
                    for obj_i = 1:size(obj.objects,1)
                        if isequal(obj.objects{obj_i,1},'U235')
                            last_r = obj.r235;
                        else
                            last_r = obj.r238;
                        end
                        [xout,~] = circcirc(s0(1),s0(2),obj.r235,obj.objects{obj_i,2},obj.objects{obj_i,3},last_r);
                        if ~isnan(xout)
                            tem_nan = true;
                            break;                
                        end
                    end
                end    
                if (~tem_nan)  
                    I_m = image(obj.u235_file,'parent', obj.Ax, 'AlphaData', obj.u235_alpha, 'XData', [s0(1)-obj.r235,s0(1)+obj.r235], 'YData', [s0(2)+obj.r235,s0(2)-obj.r235]);
                    obj.objects = [obj.objects;{'U235',s0(1),s0(2),0,0,I_m}];
                    obj.compute_obj = [obj.compute_obj;true];
                end


            end
            
            %create U238
            while size(obj.objects,1) < (obj.N235 + obj.N238)
                tem_nan = false;
                s0 = rand(1,2).*obj.gridSize; 

                if ~isempty(obj.objects)
                    for obj_i = 1:size(obj.objects,1)
                        if isequal(obj.objects{obj_i,1},'U235')
                            last_r = obj.r235; 
                        else
                            last_r = obj.r238;
                        end
                        [xout,~] = circcirc(s0(1),s0(2),obj.r238,obj.objects{obj_i,2},obj.objects{obj_i,3},last_r);
                        if ~isnan(xout)
                            tem_nan = true;
                            break;                
                        end
                    end
                end    
                if (~tem_nan)    
                    I_m = image(obj.u238_file,'parent', obj.Ax, 'AlphaData', obj.u238_alpha, 'XData', [s0(1)-obj.r238,s0(1)+obj.r238], 'YData', [s0(2)+obj.r238,s0(2)-obj.r238]);    
                    new_obj = {'U238',s0(1),s0(2),0,0,I_m};   
                    obj.objects = [obj.objects;new_obj];
                    obj.compute_obj = [obj.compute_obj;true];
                end

            end
            %create neutrons
            for i = 1:obj.Nn
               s0 = [-300,obj.gridSize(2)/2];
               v0x = obj.v0_n;%.*rand(1,2);
               p = plot(obj.Ax,s0(1),s0(2),'.k','MarkerSize',10);
               new_obj = {'n',s0(1),s0(2),v0x(1),v0x(2),p};
               obj.objects(end+1,:)=new_obj; 
               obj.compute_obj = [obj.compute_obj;true];
            end
            
            obj.Ax.XLim = [-400,1.1*obj.gridSize(1)];
            obj.Ax.YLim = [-0.1*obj.gridSize(2),1.1*obj.gridSize(2)];
            axis(obj.Ax, 'equal');
            axis(obj.Ax, 'manual');  
            set(obj.Ax,'xtick',[]);
            set(obj.Ax,'xticklabel',[]);
            set(obj.Ax,'ytick',[]);
            set(obj.Ax,'yticklabel',[]);
        end

        function obj = NewNeutronStart(obj)
            s0 = [-300,obj.gridSize(2)/2];
           v0x = obj.v0_n;%.*rand(1,2);
           p = plot(obj.Ax,s0(1),s0(2),'.k','MarkerSize',10);
           new_obj = {'n',s0(1),s0(2),v0x(1),v0x(2),p};
           obj.objects(end+1,:)=new_obj; 
           obj.compute_obj(end+1,:) = true;
        end
        
        function obj = RemoveElementsOutside(obj)
            
            objects_positions = [[obj.objects{:,2}]',[obj.objects{:,3}]'];

            XLIM = obj.Ax.XLim;
            YLIM = obj.Ax.YLim;

            obj_names = {'n','Kr91','Ba142'};
            obj_raios = [obj.rn,obj.rKr,obj.r238];

            remove_total = [];
            for obj_id = 1:3

                obj_indices = find(ismember(obj.objects(:,1),obj_names{obj_id}));
                remove_ind = or((objects_positions(obj_indices,1)+obj_raios(obj_id))<XLIM(1),...
                    or((objects_positions(obj_indices,1)-obj_raios(obj_id))>XLIM(2),...
                    or((objects_positions(obj_indices,2)+obj_raios(obj_id))<YLIM(1),...
                    (objects_positions(obj_indices,2)-obj_raios(obj_id))>YLIM(2))));

                remove_total = [remove_total;obj_indices(logical(remove_ind))];
            end
            if ~isempty(remove_total)    
                delete([obj.objects{remove_total,6}]);
                obj.objects(remove_total,:) = [];
                obj.compute_obj(remove_total,:) = [];
                
            end
        end

        function obj = UpdatePositions(obj)
            objects_table = obj.objects;
            dT = obj.dt;
            
            objToMoveInd = find(ismember(objects_table(:,1),{'n','Kr91','Ba142'}));
            wrapper = @(x,sx,sy,vx,vy) move_obj(x,sx,sy,vx,vy,dT) ;
            [images,sx,sy] = cellfun( wrapper, objects_table(objToMoveInd,6),...
                objects_table(objToMoveInd,2),objects_table(objToMoveInd,3),...
                objects_table(objToMoveInd,4),objects_table(objToMoveInd,5), 'UniformOutput', false) ;

            objects_table(objToMoveInd,2) = sx;
            objects_table(objToMoveInd,3) = sy;
            objects_table(objToMoveInd,6) = images;
            
            obj.objects = objects_table;
%           
%             
%             for i = 1:length(objToMoveInd)
%                 pi = objToMoveInd(i);
%                 objects_table{pi,6}.XData = objects_table{pi,6}.XData + objects_table{pi,4}*dT;
%                 objects_table{pi,6}.YData = objects_table{pi,6}.YData + objects_table{pi,5}*dT;
%                 objects_table{pi,2} = objects_table{pi,2} + objects_table{pi,4}*dT;
%                 objects_table{pi,3} = objects_table{pi,3} + obj.objects{pi,5}*dT;
%             end
            
        end
        
        function obj = createU239(obj,s0)
            axes(obj.Ax);
            I_m = image(obj.u239_file,'parent',obj.Ax, 'AlphaData', obj.u239_alpha, 'XData', [s0(1)-obj.r239,s0(1)+obj.r239], 'YData', [s0(2)+obj.r239,s0(2)-obj.r239]);
            obj.objects(end+1,:) = {'U239',s0(1),s0(2),0,0,I_m};
            obj.compute_obj = [obj.compute_obj;false];
            
        end

        function obj = createKrBa(obj,s0)
            axes(obj.Ax);

            %Bário
            v0xBa = obj.v0_Ba(1)*(2*randi(2)-3);
            v0yBa = obj.v0_Ba(2)*(2*randi(2)-3);

            I_m = image(obj.Kr91_file,'parent',obj.Ax, 'AlphaData', obj.Kr91_alpha, 'XData', [s0(1)-obj.rBa,s0(1)+obj.rBa], 'YData', [s0(2)+obj.rBa,s0(2)-obj.rBa]);
            obj.objects(end+1,:) = {'Ba142',s0(1),s0(2),v0xBa,v0yBa,I_m};
            obj.compute_obj = [obj.compute_obj;false];

            %Kriptônio
            v0xKr = obj.v0_Kr(1)*(2*randi(2)-3);
            v0yKr = obj.v0_Kr(2)*(2*randi(2)-3);
            I_m = image(obj.Ba142_file,'parent',obj.Ax, 'AlphaData', obj.Ba142_alpha, 'XData', [s0(1)-obj.rKr,s0(1)+obj.rKr], 'YData', [s0(2)+obj.rKr,s0(2)-obj.rKr]);

            obj.objects(end+1,:) = {'Kr91',s0(1),s0(2),v0xKr,v0yKr,I_m};
            obj.compute_obj = [obj.compute_obj;false];
            
        end

        function obj = create3Neutrons(obj,s0)
            axes(obj.Ax);
            for i = 1:3
                v0x = (2*mean(obj.v0_n))*(2*randi(2)-3);
                v0y = (2*mean(obj.v0_n))*(2*randi(2)-3);
                obj.objects(end+1,:) = {'n',s0(1),s0(2),v0x,v0y,plot(obj.Ax,s0(1),s0(2),'.k','MarkerSize',10)};
                obj.compute_obj = [obj.compute_obj;true];
            end
        end

        function [obj,TimeNextCollision,NNextCollision,UnextCollision] = nextEvent(obj)
            Allobjects = obj.objects;
            ComputeObjs = obj.compute_obj;
            
            NNextCollision = -1;
            TimeNextCollision = inf;
            UnextCollision = -1;

            index_neutrons = find(logical(ismember(Allobjects(:,1),'n').*ComputeObjs));
            index_U = find(ismember(Allobjects(:,1),{'U235','U238'}));
            for nid = 1:length(index_neutrons) 
                consider_in_future = false;
                for uid = 1:length(index_U)
                   dx = Allobjects{index_neutrons(nid),2}-Allobjects{index_U(uid),2};
                   dy = Allobjects{index_neutrons(nid),3}-Allobjects{index_U(uid),3};

                   dvx = Allobjects{index_neutrons(nid),4};
                   dvy = Allobjects{index_neutrons(nid),5};
                   dv = dvx*dvx + dvy*dvy;

                   scal = dvx*dx+dvy*dy;
                   ups = (scal*scal)- dv*(dx*dx+dy*dy-4*obj.rn*obj.r235);
                   if(ups>0)&&(scal<0)
                       deltaT = (-scal-sqrt(ups))/(dvx^2+dvy^2);
                       if ((deltaT<TimeNextCollision)&&deltaT>0)

                           NNextCollision = index_neutrons(nid);
                           TimeNextCollision = deltaT;
                           UnextCollision = index_U(uid);
                       end
                       consider_in_future = true;
                   end
                end
                ComputeObjs(nid) = consider_in_future;
            end
            obj.compute_obj = ComputeObjs;
        end
    end
end