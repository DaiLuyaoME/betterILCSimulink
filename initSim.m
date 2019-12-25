%% model parameters
modelTypeName = {'rigidBody','doubleMassNonColocated','doubleMassColocated'};
modelInfo.mass = [5,20];
modelInfo.fr = 200;
modelInfo.dampRatio = 0.03;
modelInfo.type = modelTypeName{2};
fs = 5000;
Ts = 1/fs;
Gp = createPlantModel(modelInfo);

%% delay factor
delayCount = 1.3;
s = tf('s');
delayModel = exp(-delayCount*Ts*s);
delayModel = pade(delayModel,2);

%% generate plant model with delay
GpWithDelay = Gp * delayModel;
GpDis = c2d(GpWithDelay,Ts,'zoh');


%% ideal feedforward coefficients 
m = modelInfo.mass;
% idealAccCoef = sum(m);
% idealJerkCoef = sum(m) * tau;
% idealSnapCoef = sum(m) * ( 1/wn.^2 + 0.5 * tau.^2);
%%
sigma = 2;%�����ı�׼���λm
varNoise=sigma*sigma;%ע�⣬��������ģ���е�Noise Power ��Ҫ���varNoise*Ts
noisePower=varNoise*Ts;
%%
A1 = 8;
A2 = 5;
A3 = 0;