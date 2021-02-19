%%
% This is the static penalty case
% To change into dynamic penalty case, uncomment the line 98.
% To get the best results for each run, open the 'fff' variable
% To get the corresponding best variables, open the 'rgbest' variable
% To get best result for each iteration and each run, open 'ffmin' variable

% The explanation for the code is given in a section wise below

%%
clear all; close all; clc;
warning off;
tic

%% Input Parameters
% Varies for each cases 
lb = [0 0];            % lower bound for nv = 2
ub = [50 50];             % upper bound for nv = 2

ns = 60;                % population of swarm particles
nv = 2;                 % number of design variables
% change nv if number of varible changes

iter = 500;             % iterations total
func = @P5;             % assigning variable to call to the function
                        
maxrun = 10;            % total number of runs
c1 = 1.8;               % acceleration parameter (local)
c2 = 2.2;               % acceleration parameter (global)
w=1.3;                  % inertia weight
delt = 1;               % delta t (time step) value
pfc = 1.1;              % penalty function value

%% Starting the PSO for all the runs
% initiating the 'for loop' for obtaining the results for all the 10 runs 

for run = 1:maxrun
    
    % initializing the position and velocity of the swarm
    x = zeros(ns,nv);           
    random_n = rand(ns,nv);     % calculates the random number matrix
    for i = 1:nv
        x(:,i) = lb(1,i) + (ub(1,i)-lb(1,i)).*random_n(:,i);
    end
    v = zeros(ns,nv);
   
    %% calculating the first iteration to get the local and global best
    % position of all the particles
    % ValF0 indicates the fitness function value
    % Valbest indicates the best fitness function value
    % xn0 indicates the new position matrix of swarm
    % xbest is the best local position of each particle
    % gbest is the best global position among all the particles
    
    for i = 1:ns
        [ValF0(i,1)] = func(x(i,:),pfc);
        xn0(i,:) = x(i,:);
    end
    [Valbest,index0] = min(ValF0);
    xbest = xn0;
    gbest = xn0(index0,:);
    
    it = 1;     % iteration count initiated to 1
    w = 1.3;    % after every 1 run, the inertia weight reassign to its original value
    pfc = 1.1;  % after every 1 run, the penalty function value reassign to its original value
    
    %% running the main PSO loop for 'iter' iterations
    
    % first updating the particle velocity using the given equation
    while it<=iter
        for i = 1:ns
            for j = 1:nv
                r1 = rand();
                r2 = rand();
                v(i,j) = w*v(i,j) + c1*r1*((xbest(i,j)-x(i,j))/delt) + (c2*r2*(gbest(1,j)-x(i,j))/delt);
            end
        end
        
        % now updating the particle best position 
        for i = 1:ns
            for j = 1:nv
                x(i,j) = x(i,j) + v(i,j)*delt;
            end
        end
        
        % ensuring that the particle does not go out of the design space
        for i = 1:ns
            for j = 1:nv
                if x(i,j) > ub(1,j)
                    x(i,j) = ub(1,j);
                end
                if x(i,j) < lb(1,j)
                    x(i,j) = lb(1,j);
                end
            end
        end
        
%         pfc = pfc*1.1;          % dynamic penalty function
        
        % calculating the fitness function value
        % when calling the function, the input penalty function 'pfc' is
        % used for constraint one only
        for i = 1:ns
            [ValF(i,1)] = func(x(i,:),pfc);
        end
        
        % comparing the fitness function value with its previous best value
        % and updating the same if new value is less than the old and 
        % also it updates the best position xbest
        
        for i = 1:ns
            if ValF(i,1) < ValF0(i,1)
                xbest(i,:) = x(i,:);
                ValF0(i,1) = ValF(i,1);
            end
        end
        
        % finding and storing the best value for each run and iterations
        % fmin indicates the best particle fitness value for the running
        % particular iteration and sorts out from swarm size
        % ffmin stores the best fitness value for each run and iteration
        % ffite shows the iteration performed for each run - eventually it
        % will be the 'maxrun' x 1 matrix with each element equal to 'iter'
        
        [fmin,index]=min(ValF0);       % finding out the best particle
        ffmin(it,run)=fmin;            % storing best fitness
        ffite(run)=it;                 % storing iteration count
       
        % updating gbest and best fitness function value
        if fmin<Valbest
            gbest=xbest(index,:);
            Valbest=fmin;
        end
        
        %% displaying iteration results initiation
        
        if it==1
            fprintf('Iteration    Best particle    Objective fun\n');
        end
        fprintf('%8g  %8g          %8.4f\n',it,index,Valbest);
        
        it=it+1;        % updating the iteration count 
        w = w*0.999;    % damping factor for inertia weight
    end
    
    %% Plotting the convergence characteristics for each run
    % Also saves the plot with the given name in temp
    
    figure;
    plot(ffmin(1:ffite(run),run),'-k');
    temp=['P5_c1_',num2str(c1),'_c2_',num2str(c2),'_',num2str(run),'.png'];
    saveas(gcf,temp);
    
    % PSO main program ------------------------------------ends
   
    %% Calculating the original function value by using the gbest (global best position) for each run
    % fvalue indicates the best function value
    % fff stores the best function value for each run
    % rgbest stores the best function variable for each run
    
    m = 1;
    k = gbest(1);
    c = gbest(2);
    omega = 0.1;
    F0 = 1;
    
    C = sqrt( (k-m*omega^2)^2 + (c*omega)^2 );
    alpha = atan( c*omega/(k -m*omega^2));
    omegad = sqrt( k/m - (c/(2*m))^2 );
    A = - (F0/C) * cos(alpha);
    B = -(F0/(C*omegad)) * ( omega*sin(alpha) + c/(2*m) * cos(alpha));
    
    f = fopen('MeasuredResponse.dat');
    data1 = textscan(f, '%f %f');
    data1 = cell2mat(data1); % convert to matrix from cell array
    fclose(f);
    
    t = 0;
    for i = 1:500
        term1 = A*cos(omegad*t) + B*sin(omegad*t);
        term2 = exp( - c*t/(2*m) );
        term3 = (F0/C)*cos(omega*t - alpha);
        u(i) = term1*term2 + term3;
        time(i) = t;
        t = t + 10/499;
        udot(i) = term2*(B*omegad*cos(omegad*t) - A*omegad*sin(omegad*t)) + (A*cos(omegad*t) + B*sin(omegad*t))*term2*(-c/(2*m)) - (F0/C)*omega*sin(omega*t-alpha);
        udotdot(i) = -term2*omegad^2*term1 - (c/(2*m))*term2*B*omegad*cos(omegad*t) + term2*(c/(2*m))*A*omegad*sin(omegad*t) + ((c/(2*m))^2)*term1*term2 - (c/(2*m))*term2*omegad*(B*cos(omegad*t)-A*sin(omegad*t)) - (F0/C)*(omega^2)*cos(omega*t-alpha);
        s1(i) = m*udotdot(i) + c*udot(i) + k*u(i) - F0*cos(omega*t);
        ans1(i) = (abs((u(i)'-data1(i,2))));
        ans1(i) = ans1(i).^2;
    end
    fvalue = sum(ans1);
    
    fff(run)=fvalue;
    rgbest(run,:)=gbest;
    fprintf('--------------------------------------\n');
end
toc
% optimization -------------------ends for all the runs

%% Displaying the final best result
% displays the best function value and variables among all the runs

fprintf('\n\n');
fprintf('*****************************************************\n');
fprintf('Final Results-------------------------\n');
[bestfun,bestrun]=min(fff)
best_variables=rgbest(bestrun,:)
fprintf('*********************************************************\n');

% displaying the convergence plot for the best run

plot(ffmin(1:ffite(bestrun),bestrun),'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('PSO convergence characteristic')
saveas(gcf,'Fit_BestP5.png');

%##############################################-----------------end