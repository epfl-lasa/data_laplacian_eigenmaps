clear all;
close all;

%% Saving the generated data
% Structure of the data = [X alpha Xdot alpha_dot Time Subdynamic_Id]
% By default, subdynamic_Id is 1. Changes in case there are multiple
% subdynamics in the demonstrations

% Params for the demonstrations
demoCount = 1; % (Data's Dimensionality + 1)
% c = [1, 3, 7, 10, 15, 20]; % Swirl
c = [1]; % Swirl

% N = [100,500,1000,2000]; % Point count per demonstration
N = [500]; % Point count per demonstration

% Other params
label = 1;

k = 1.4; % Speed of spiral traversal
thetaDot = 2;
dT = 0.006;
attractor = [0 0];

pointCount = 10;
initPoints = zeros(pointCount, 2);
initPoints(:,1) = -3 + (3+3)*rand(pointCount,1);
initPoints(:,2) = -3 + (3+3)*rand(pointCount,1);



for nIndex = 1:length(N)
  for cIndex = 1:length(c)
    demo_struct = {'position', 'velocity', 'time', 'labels'};
    demo = cell(pointCount, 1);
    nVal = N(nIndex);
    cVal = c(cIndex);
    
    %% Generating the data
    for count = 1:pointCount
      startPoint = initPoints(count, :);
      genData = [startPoint];
      currentPoint = startPoint;
      
      while norm(currentPoint - attractor) > 1e-5

        tempVel = spiralVel(currentPoint, k, thetaDot);
        currentPoint = currentPoint + tempVel*dT; 
%         genData = [genData; [tempPos, tempVel, t]];
        genData = [genData; [currentPoint]];
      end    
      demo{count} = genData;
    end
    
    
%     name = strcat("N_",int2str(nVal),"_c_",int2str(cVal),".mat");
%     save(name ,'demo','demo_struct', "-mat");
%     
    for count = 1:pointCount
        demo{count}(1,:)
        plot(demo{count}(:,1), demo{count}(:,2), 'r')
        hold on
    end
  end
end


%% Helper functions

function vel = spiralVel(currentPos, k, thetaDot)
 
  x = currentPos(1);
  y = currentPos(2);
  
  xDot = -k*x - thetaDot*y;
  yDot = -k*y + thetaDot*x;

  vel = [xDot yDot];
end










