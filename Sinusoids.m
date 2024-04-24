clear all;
close all;

%% Helper functions

%% Saving the generated data
% Structure of the data = [X Xdot Time Subdynamic_Id]
% By default, subdynamic_Id is 1. Changes in case there are multiple
% subdynamics in the demonstrations

% Params for the demonstrations
demoCount = 1; % (Data's Dimensionality + 1)
N = 500; % Point count per demonstration
% Other params
label = 1;
cList = [1, 3, 5, 7, 10];
K = 5;
T = 1.5;
t = linspace(0,T,N);
dim = 10;

addNoise = false;

rng(0);
A = rand(1, dim); % Amplitude

rng(1);
phi = pi/2*rand(1, dim); % Initial Phase

rng(2);
f = rand(1, dim);
freq = {}; % Frequency for difference complexities
for cIndex = 1:length(cList)
    c = cList(cIndex);
    freq{cIndex} = (c - 1)*f + 1;
end

for cIndex = 1:length(cList)
    genData = [];
    for n = 1:N
        demo_struct = {'position', 'velocity', 'time', 'labels'};
        demo = cell(demoCount, 1);

        % Generating the data
        [theta, thetaDot] = getTheta(K, t(n));
        tempPos = A.*sin(phi + freq{cIndex}*theta);
        if addNoise
            tempPos = awgn(tempPos, 25, 'measured');
        end
        tempVel = thetaDot*A.*freq{cIndex}.*cos(phi + freq{cIndex}*theta);
        genData = [genData; [tempPos, tempVel, t(n)]];
    end    
    genData = [genData ones(size(genData,1),1)];
    demo{1} = genData';

    name = strcat("N_",int2str(N),"_c_",int2str(cList(cIndex)),".mat");
    save(name ,'demo','demo_struct', "-mat");
end
%% Visualising the data

figure;
for demoIndex = 1:demoCount
   xlabel('t')
   ylabel('y')
   for d = 1:dim
       hold on
       plot(t, demo{1}(d,:))
   end
   hold off
end

%% Helper functions

function [theta, thetaDot] = getTheta(k, t)
    theta = pi/2*(1 - exp(-k*t)); 
    thetaDot = k*pi/2*exp(-k*t);
end
