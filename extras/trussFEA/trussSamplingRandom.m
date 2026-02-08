function [relV,relJ] = trussSamplingRandom(t,nSamples)
% Relative volume and relative compliance for trusses with random areas
N = t.myNumTrussBars;
A0 = t.myArea;
relV = zeros(1,nSamples);
relJ = zeros(1,nSamples);
V0 = t.getVolume();
J0 = t.getCompliance();
AMinRel = 0.1;
AMaxRel = 2;
for i = 1:nSamples   
    A = A0.*(AMinRel + (AMaxRel-AMinRel)*rand(1,N)); %
    t = t.assignA(A);
    t = t.assemble();
    t = t.solve();
    relV(i) = t.getVolume()/V0;
    relJ(i) = t.getCompliance()/J0;
end
figure; plot(relV,relJ,'*'); hold on;
xlabel('Relative Volume');
ylabel('Relative Compliance');
grid on;
plot(1,1,'dr')

