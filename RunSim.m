%%
% accCoef = accCoefLast;
% jerkCoef = jerkCoefLast;
% snapCoef = snapCoefLast;

accCoef = 25;
jerkCoef = 0.0090;
snapCoef = 2.415e-6 + 4.5094e-07;

% accCoef = 25;
% jerkCoef = 0;
% snapCoef = 0;



trajParameters.dis = 0.04;
trajParameters.vel = 0.25;
trajParameters.acc = 10; 
trajParameters.jerk = 800;
trajParameters.snap = 64000;

alpha =  0;

sim('main',[0 0.02]);