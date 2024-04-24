clear
close all;

%% Saving the generated data
% Structure of the data = [X alpha Xdot alpha_dot Time Subdynamic_Id]
% By default, subdynamic_Id is 1. Changes in case there are multiple
% subdynamics in the demonstrations

generateRandom = true; % Generate trajectories from random starting points or not
saveData = false; % Save data as .mat

% Swirl
% cList = [1, 2, 3, 5, 7, 10];
cList = [3];


% Point-Count
NList = [500];

% Speed of spiral traversal
k1 = 0.03; 
k2 = 0.03;

% dT = 0.002;
dT = 0.003;

limit = 100000; % Max number of points to be generated (worst case)
attractor = [0 0 -1];

%% Random starting point computation
pointCount = 10;
initPoints = zeros(pointCount, 3);
initPoints(:,1) = -0.3 + (0.3+0.3)*rand(pointCount,1);
initPoints(:,2) = -0.3 + (0.3+0.3)*rand(pointCount,1);
initPoints(:,3) = sqrt(1 - (initPoints(:,1).^2 + initPoints(:,2).^2));
%%
initPoints
%% Generating the data
% Generating forward orbits for visualisation
for cIndex = 1:length(cList)
    for nIndex = 1:length(NList)
        N = NList(nIndex);
        c = cList(cIndex);
        disp([">> c Value: ", c, " N Value: ", N])
        if generateRandom
            demoR = cell(pointCount, 1);
            for count = 1:pointCount

              startPoint = initPoints(count, :);
              genData = [startPoint];
              currentPoint = startPoint;

              psiOld = atan2(startPoint(2),startPoint(1));
              turn = 0;
              if psiOld < 0
                  psiOld = psiOld + 2*pi;
              end


              while norm(currentPoint - attractor) > 1e-2 && size(genData, 1) < limit
  
                [tempVel, psiOld, turn] = spiralVel(currentPoint, k1, k2, c, psiOld, turn);
                currentPoint = currentPoint + tempVel*dT; 
                genData = [genData; [currentPoint]];
              end    
              demoR{count} = genData;
            end

            for count = 1:pointCount
            %     scatter3(demoR{count}(1:10:end,1), demoR{count}(1:10:end,2), demoR{count}(1:10:end,3), 'r')
                plot3(demoR{count}(:,1), demoR{count}(:,2), demoR{count}(:,3), 'r')
                hold on
            end
        end

        % Generating and saving the demonstration
        demo_struct = {'position', 'velocity', 'time', 'labels'};
        demo = cell(1, 1);

        startPoint = [0 0 1];
        genData = [];
        currentPoint = startPoint;
        prevPoint = [10 10 10];

        psiOld = atan2(startPoint(2),startPoint(1));
        turn = 0;

        t = 0;
        while norm(currentPoint - attractor) > 1e-3 && size(genData, 1) < limit && norm(currentPoint - prevPoint) > 1e-5

            [tempVel, psiOld, turn] = spiralVel(currentPoint, k1, k2, c, psiOld, turn);
            genData = [genData; [currentPoint, tempVel, t]];
            prevPoint = currentPoint;
            currentPoint = currentPoint + tempVel*dT;
            t = t+dT;

        end
        genData(end,:)
        size(genData,1)

        
        figure;
        genData = [genData ones(size(genData,1),1)];
        plot3(genData(:,1), genData(:,2), genData(:,3), 'r')
        hold on

        % Subsample N points from the generated data
        genData  = genData (1:ceil(size(genData,1)/N):end, :);
        scatter3(genData(:,1), genData(:,2), genData(:,3), 'b')
        hold off

        % Saving the data
        demo{1} = genData';
        if saveData
            name = strcat("N_",int2str(N),"_c_",int2str(c),".mat");
            save(name ,'demo','demo_struct', "-mat");
        end
    end
end
%% Helper functions

function [vel, psiN, turnN] = spiralVel(currentPos, k1, k2, c, psiOld, turn)
 
  x = currentPos(1);
  y = currentPos(2);
  z = currentPos(3);
  
  theta = acos(z);
  psi = atan2(y,x);
  if psi < 0
      psi = psi + 2*pi;
  end

  psiTemp = turn*2*pi + psi;
  if psiTemp < psiOld
      turn = turn + 1;
  end
  
  psiOld = psiTemp;
  psi = turn*2*pi + psi;

  thetaDot = k1*(pi - theta);
  psiDot = k2*(2*c*pi - psi);
  
  xDot = -sin(psi)*sin(theta)*psiDot + cos(psi)*cos(theta)*thetaDot;
  yDot = cos(psi)*sin(theta)*psiDot + sin(psi)*cos(theta)*thetaDot;
  if z == 1
      zDot = -sqrt(x^2+y^2)*thetaDot;
  else
      zDot = -sin(theta)*thetaDot;
  end


  vel = [xDot yDot zDot];
  psiN = psiOld;
  turnN = turn;
end









