close all;clear;
gridSize = [800,800];
N235 = 100;
N238 = 0;
Nn = 1; 
r235 = 15;
r238 = 15;
r239 = 15;
rn = 5;
rBa = 10;
rKr = 10;

v0_n = [1000,0];
v0_Kr = [150,150];
v0_Ba = [150,150];

t0 = 0;
tf = 30;
dt = 0.0001;

%keep track:
%iniciar modelo:   
%model
gamma_param = 4;
u_model = N235;
n_model = Nn;
m_model = N238;
t_model = t0;
for i = t0+dt:dt:0.5
    dn = (3*n_model(end)*u_model(end)-gamma_param*n_model(end)*m_model(end))*dt;
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


figure;
plot(t_model(t_model<0.02),n_model(t_model<0.02),'r','LineWidth',3);
ylabel('quantidade de nêutrons')
xlabel('tempo')


plot(t_model,u_model,'k','LineWidth',3);
hold on
plot(t_model(t_model<0.02),n_model(t_model<0.02),'r','LineWidth',3);
%plot(t_model,m_model,'g');
hold off;
%legend({'u';'n';'m'})



Fig = figure('units','normalized','outerposition',[0 0 1 1]);
Ax = gca;
Chain = ChainReaction(Ax);
Chain = Chain.GetParameters(gridSize, N235, N238, Nn, r235, r238, r239, rn, rBa, rKr, v0_n, v0_Kr, v0_Ba);
Chain = Chain.GenerateElements();
Chain = Chain.ChainSim(t0,tf,dt);

8.05*3120+9*72+10*72+7*108+10*72+9*108;
