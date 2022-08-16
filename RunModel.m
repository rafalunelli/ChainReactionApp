function [n_model,u_model,m_model,t_model] = RunModel(N235,Nn,N238)
gamma_param = 4;

u_model = N235;
n_model = Nn;
m_model = N238;
t_model = 0;

dt = 0.0001;

for i = dt:dt:0.5
    dn = (2*n_model(end)*u_model(end)-gamma_param*n_model(end)*m_model(end))*dt;
    du = (-n_model(end)*u_model(end))*dt;
    dm = (-gamma_param*n_model(end)*m_model(end))*dt;
    u_model(end+1) = u_model(end)+du;
    n_model(end+1) = n_model(end)+dn;
    m_model(end+1) = m_model(end)+dm;
    t_model(end+1) = t_model(end)+dt;
end

idx = logical((u_model>0.001).*(n_model<(3*N235+Nn)));

t_model = t_model(idx);
m_model = m_model(idx);
n_model = n_model(idx);
u_model = u_model(idx);

end