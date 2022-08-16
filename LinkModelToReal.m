function [n_model,u_model,m_model,t_model] = LinkModelToReal(v0_n,n_real,u_real,m_real,t_real,n_model,u_model,m_model,t_model)


%tempo
t_model = t_model - min(t_model);
t_model = (max(t_real)-min(t_real))*t_model/max(t_model) + min(t_real);

% 
% %uranio 235
% u_model = u_model - min(u_model);
% u_model = (max(u_real)-min(u_real))*u_model/max(u_model) + min(u_real);
% 
% %uranio 238
% m_model = m_model - min(m_model);
% m_model = (max(m_real)-min(m_real))*m_model/max(m_model) + min(m_real);
% 
% %uranio 238
% n_model = n_model - min(n_model);
% n_model = (max(n_real)-min(n_real))*n_model/max(n_model) + min(n_real);
%             
end