%% fit the model1
load parameters.mat
load g0c.mat
load g1c.mat
load g2c.mat

model0 = fitrnet(para(1:800,:),g0c(1:800),"Standardize",true,"LayerSizes",[20 20 20], ...
    'IterationLimit',2000,'LossTolerance',1e-12);
model1 = fitrnet(para(1:800,:),g1c(1:800),"Standardize",true,"LayerSizes",[20 20 20], ...
    'IterationLimit',2000,'LossTolerance',1e-12);
model2 = fitrnet(para(1:800,:),g2c(1:800),"Standardize",true,"LayerSizes",[20 20 20], ...
    'IterationLimit',2000,'LossTolerance',1e-12);

save('g0_ml.mat','model0')
save('g1_ml.mat','model1')
save('g2_ml.mat','model2')

loss(model1,para(801:1000,:),g1c(801:1000))

%% plotting
figure('Position', [100, 100, 1000, 400])
FS = 'fontsize'; FW = 'fontweight'; NO = 'normal'; LW = 'linewidth';
subplot(2,1,1)
plot(g1c(1:800),'kx','DisplayName','LARS'); hold on 
% plot(para(1:800,:)*coeff,'r-o','DisplayName','Neural network training');
plot(predict(model1,para(1:800,:)),'ro','DisplayName','Neural network training');
xlabel('Number', FS,14)
ylabel('$G_1$','interpreter','latex', FS,14)
title('Model training', FS,14)
legend
set(gca,'PlotBoxAspectRatio',[7 1 1])
ylim([0 .08])
subplot(2,1,2)
plot(g1c(801:1000),'kx','DisplayName','LARS'); hold on 
% plot(para(801:1000,:)*coeff,'b-o','DisplayName','Neural network prediction');
plot(predict(model1,para(801:1000,:)),'bo','DisplayName','Neural network prediction');
xlabel('Number', FS,14)
ylabel('$G_1$','interpreter','latex', FS,14)
title('Model validation', FS,14)
legend
set(gca,'PlotBoxAspectRatio',[7 1 1])
ylim([0 .08])
exportgraphics(gcf,'ml_smooth.pdf','Resolution',600)