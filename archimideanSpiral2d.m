clear all;
close all;

%% Saving the generated data
% Structure of the data = [X Xdot Time Subdynamic_Id]
% By default, subdynamic_Id is 1. Changes in case there are multiple
% subdynamics in the demonstrations

% Params for the demonstrations
demoCount = 1;

c = [2]; % Swirls in spiral
r = 0.1; % Radius of sphere
T = 12; % Time to make the spiral
N = [500]; % Point count per demonstration

% Other params
label = 1;

for nIndex = 1:length(N)
  for cIndex = 1:length(c)
    demo_struct = {'position', 'velocity', 'time', 'labels'};
    demo = cell(demoCount, 1);
    nVal = N(nIndex);
    cVal = c(cIndex);
    
    target = cVal*2*pi;
    
    dT = T/nVal; % Recording frequency
    
    %% Generating the data
    for count = 1:demoCount
      t = 0;
      genData = [];
      for pointCount = 1:nVal
        tempPos = spiralPos(t,r,cVal,T);
        tempVel = spiralVel(t,r,cVal,T);
        genData = [genData; [tempPos, tempVel, t]];
        t = t + dT;
      end    
      genData = [genData ones(size(genData,1),1)];
      demo{count} = genData';
    end
    
    name = strcat("N_",int2str(nVal),"_c_",int2str(cVal),".mat");
    save(name ,'demo','demo_struct', "-mat");
      
  end
end
%% Visualising the data

figure
for demoIndex = 1:demoCount
   plot(demo{demoIndex}(1,:), demo{demoIndex}(2,:), 'r')
   xlabel('x')
   ylabel('y')
end



%% Helper functions
function pos = spiralPos(t,r,c,T)
    target = 2*pi*c;
    thetaDot = target/T;
    theta = 0 + t*thetaDot;
    % rad = b * theta, b = R/target
    b = r/target;
    rad = 0 + b*theta;
    
    x = rad*cos(theta);
    y = rad*sin(theta);
    
    pos = [x y];
end

function vel = spiralVel(t,r,c,T)
    target = 2*pi*c;
    thetaDot = target/T;
    theta = 0 + t*thetaDot;
    b = r/target;
    
    xDot = b*thetaDot*cos(theta) - b*theta*thetaDot*sin(theta);
    yDot = b*thetaDot*sin(theta) + b*theta*thetaDot*cos(theta);

    vel = [xDot yDot];
end

%%