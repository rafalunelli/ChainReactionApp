classdef ChainPlotClass
    
    properties
        ax_real
        ax_model
        
        n_real
        u_real
        m_real
        t_real
        
        n_model
        u_model
        m_model
        t_model
        
    end
    
    methods
        function obj = ChainPlotClass(ax_real,ax_model)
            
            obj.ax_real = ax_real;
            obj.ax_model = ax_model;
        end
        
        function obj = getData_real(obj,n_real,u_real,m_real,t_real)
            obj.n_real = n_real;
            obj.u_real = u_real;
            obj.m_real = m_real;
            obj.t_real = t_real;
        end
        
        function obj = getData_model(obj,n_model,u_model,m_model,t_model)
            obj.n_model = n_model;
            obj.u_model = u_model;
            obj.m_model = m_model;
            obj.t_model = t_model;
        end
        
        function plotData(obj,plot_mode)
            
            switch plot_mode
                case 'Nêutrons'
                    plot(obj.ax_real,obj.t_real,obj.n_real);
                    obj.ax_real.YLabel.String = 'Quantidade: nêutrons';
                    obj.ax_model.YLabel.String = 'Quantidade: nêutrons';
                    plot(obj.ax_model,obj.t_model,obj.n_model);
                case 'Urânio 235'
                    plot(obj.ax_real,obj.t_real,obj.u_real);
                    obj.ax_real.YLabel.String = 'Quantidade: U-235';
                    obj.ax_model.YLabel.String = 'Quantidade: U-235';
                    plot(obj.ax_model,obj.t_model,obj.u_model);
                case 'Urânio 238'
                    plot(obj.ax_real,obj.t_real,obj.m_real);
                    obj.ax_real.YLabel.String = 'Quantidade: U-238';
                    obj.ax_model.YLabel.String = 'Quantidade: U-238';
                    plot(obj.ax_model,obj.t_model,obj.m_model);
            end
            obj.ax_real.XLabel.String = 'Tempo';
            obj.ax_model.XLabel.String = 'Tempo';
            
            obj.ax_real.Title.String = 'Resultado: Simulação';
            obj.ax_model.Title.String = 'Resultado: Modelização';
        end
    end
end

