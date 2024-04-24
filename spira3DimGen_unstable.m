clear all;
close all;

%% Saving the generated data
% Structure of the data = [X Xdot Time Subdynamic_Id]
% By default, subdynamic_Id is 1. Changes in case there are multiple
% subdynamics in the demonstrations

% Params for the demonstrations
demoCount = 1;
% c = [1,3,7,10,15]; % Swirl
% r = 1; % Radius of sphere
% k = 1; % Speed of spiral traversal
% N = [500]; % Point count per demonstration

% c = [2]; % Swirl
% r = 0.15; % Radius of sphere
% k = 0.1; % Speed of spiral traversal
% N = [500]; % Point count per demonstration

c = [4]; % Swirl
r = 0.2; % Radius of sphere
k = 0.1; % Speed of spiral traversal
N = [500]; % Point count per demonstration

% Other params
label = 1;

for nIndex = 1:length(N)
  for cIndex = 1:length(c)
    demo_struct = {'position', 'velocity', 'time', 'labels'};
    demo = cell(demoCount, 1);
    nVal = N(nIndex);
    cVal = c(cIndex);
    dT = pi/(nVal*k); % Recording frequency
    
    %% Generating the data
    for count = 1:demoCount
      t = 0;
      genData = [];
      for pointCount = 1:nVal
        tempPos = spiralPos(t,r,cVal,k);
        tempVel = spiralVel(t,r,cVal,k);
        genData = [genData; [tempPos, tempVel, t]];
        t = t + dT;
      end    
      genData = [genData ones(size(genData,1),1)];
      demo{count} = genData';
    end
    
    name = strcat("N_",int2str(nVal),"_c_",int2str(cVal),".mat");
    save(name ,'demo','demo_struct', "-mat");
   
   
   
    demoIndex = 1;
    subplot(length(c),length(N),(cIndex-1)*length(N)+nIndex);
    plot3(demo{demoIndex}(1,:), demo{demoIndex}(2,:), demo{demoIndex}(3,:), 'r')
    title(strcat("N: ", int2str(nVal), " C: ", int2str(cVal)));
  end
end
%% Visualising the data

% figure
% for demoIndex = 1:demoCount
%    plot3(demo{demoIndex}(1,:), demo{demoIndex}(2,:), demo{demoIndex}(3,:), 'r')
%    xlabel('x')
%    ylabel('y')
%    zlabel('z')
%    hold on
% end



%% Helper functions
function pos = spiralPos(t,r,c,k)
  T = k*t;
  x = r*sin(T)*cos(c*T);
  y = r*sin(T)*sin(c*T);
  z = r*cos(T);
  pos = [x y z];
end

function vel = spiralVel(t,r,c,k)
  T = k*t;
  xDot = r*cos(T)*cos(c*T)-r*c*sin(T)*sin(c*T);
  yDot = r*cos(T)*sin(c*T)+r*c*sin(T)*cos(c*T);
  zDot = -r*sin(T);
  vel = [xDot yDot zDot];
end

%%