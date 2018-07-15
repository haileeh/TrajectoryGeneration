%% plot accelerations
load('DV_struct_compare.mat');
m = 5;
DV_array.Tend = DV_computed_struct.Tend;
DV_array.Alin = DV_computed_struct.Alin;
DV_array.Acen = DV_computed_struct.Acen;
DV_array.Acor = DV_computed_struct.Acor;
DV_array.Aang = DV_computed_struct.Aang;
DV_array.DV_out = DV_computed_struct.DV_out;

%% sum of norms:
ALIN_array = zeros(100,5);
ACEN_array = ALIN_array;
ACOR_array = ALIN_array;
AANG_array = ALIN_array;

h = length(DV_array.Tend);
for j=1:h
    for i=1:100
        ALIN_array(i,j) = norm(DV_array.Alin(:,i,j));
        ACEN_array(i,j) = norm(DV_array.Acen(:,i,j));
        ACOR_array(i,j) = norm(DV_array.Acor(:,i,j));
        AANG_array(i,j) = norm(DV_array.Aang(:,i,j));
    end
end

ATOT = ALIN_array + ACEN_array + ACOR_array + AANG_array;
ALIN_frac = ALIN_array ./ ATOT;
ACEN_frac = ACEN_array ./ ATOT;
ACOR_frac = ACOR_array ./ ATOT;
AANG_frac = AANG_array ./ ATOT;

ARAD_frac = ALIN_frac + ACEN_frac;
ATAN_frac = ACOR_frac + AANG_frac;

for k=1:h
    t = linspace(0,DV_array.Tend(k),length(ARAD_frac));
    figure
    plot(t,ARAD_frac(:,k)); hold on; grid on;
    plot(t,ATAN_frac(:,k));
    legend('Radial','Tangential');
    str = ['SUM OF NORMS: \Delta V = ', num2str(DV_array.DV_out(k))];
    title(str);
end

%% norm of sums
arad_array = zeros(3,100,5);
atan_array = arad_array;
ARAD = zeros(100,5);
ATAN = ARAD;
ARAD_ratio = zeros(100,5);
ATAN_ratio = ARAD_ratio;

for j=1:h
    for i=1:100
        arad_array(:,i,j) = DV_array.Alin(:,i,j) + DV_array.Acen(:,i,j);
        atan_array(:,i,j) = DV_array.Acor(:,i,j) + DV_array.Aang(:,i,j);
        ARAD(i,j) = norm(arad_array(:,i,j));
        ATAN(i,j) = norm(atan_array(:,i,j));
    end
end

ARAD_ratio = ARAD ./ ATOT;
ATAN_ratio = ATAN ./ ATOT;

for k=1:h
    t = linspace(0,DV_array.Tend(k),length(ARAD_ratio));
    figure
    plot(t,ARAD_ratio(:,k)); hold on; grid on;
    plot(t,ATAN_ratio(:,k));
    legend('Radial','Tangential');
    str = ['NORM OF SUM: \Delta V = ', num2str(DV_array.DV_out(k))];
    title(str);
end

%% plot all DVs
figure
plot(DV_array.Tend, DV_array.DV_out,'or'); hold on; grid on;
xlabel('Final Time [s]'); ylabel('DV');
title('DV');

%% Vector comparison

for k=1:h
    t = linspace(0,DV_array.Tend(k),length(arad_array));
    figure();
    subplot(3,1,1)
    plot(t(1:end),arad_array(1,:,k),'b'); hold on;
    plot(t(1:end),atan_array(1,:,k),'r');
    legend('Linear1+Centripetal1','Coriolis1+Angular1');
    str = ['Vector Elements: \Delta V = ', num2str(DV_array.DV_out(k))];
    title(str);
    xlabel('t [s]')
    ylabel('Acceleration Norm')
    grid on;
    subplot(3,1,2)
    plot(t(1:end),arad_array(2,:,k),'b'); hold on;
    plot(t(1:end),atan_array(2,:,k),'r');
    legend('Linear2+Centripetal2','Coriolis2+Angular2');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
    subplot(3,1,3)
    plot(t(1:end),arad_array(3,:,k),'b'); hold on;
    plot(t(1:end),atan_array(3,:,k),'r');
    legend('Linear3+Centripetal3','Coriolis3+Angular3');
    xlabel('t [s]')
    ylabel('Acceleration Component [m/s^2]')
    grid on;
end