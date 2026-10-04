clear; close;
tic;    
hold;
    n_sim = 1000; %Number of simulation
    r = 4; L = 20; %Channel radius and height in um
    xlim([-r, r]); ylim([-L, 0]);
    pbaspect([2*r, L, 1]);
    BP = 0.36; %Binding probability
    FS = 132.7; %Average flow speed in um/s
    DC = 8e-11; %Diffusion coefficient in m2/s
    time_step = 1e-5; %Simulation time step in s
    
    MFS = 2*FS*time_step; %Convection diplacement in um for the time step at maximum flowspeed
    DL = sqrt(4*DC*time_step)*1e6; %Diffusion length in um for the time step
    captured = 0; %Capture counter
    SPTs = cell(n_sim, 1); %SPT storage
    all_particle_no_step = zeros(n_sim, 1);
    all_particle_no_collision = zeros(n_sim, 1);
    all_particle_final_point = zeros(n_sim, 2);
    
    parfor m = 1:n_sim
        particle_trace = [2*r*rand-r,0]; %Random entry of a single molecule
        location_this_step = [0,0]; %Assign an empty vector for location of step i
        i = 0; n = 0; %assign indices. i: step counter, n: collision counter
        
        while 1
            i = i + 1;
            drift = [0,(particle_trace(end,1).^2/r^2-1).* MFS];
            theta = 2 * pi * rand; % random direction drawn
            Random_walk_direction = [cos(theta), sin(theta)]; %Vectorize the direction
            DLR = abs(randn*DL); %Pick a random walk size
            location_this_step = particle_trace(i,:) + DLR .* Random_walk_direction + drift;
            location_this_step(1) = clip(location_this_step(1),-r,r);
            
            particle_trace(i+1,:) = location_this_step;
            
            if particle_trace(end,2)>0
                %Restart...
                particle_trace = [2*r*rand-r,0];
                location_this_step = [0,0];
                i = 0; n = 0;
                continue;
            elseif particle_trace(end,2) < -L
                break;
            end
            
            if or(particle_trace(end,1) == -r,particle_trace(end,1) == r)
                n = n + 1; %Counting number of collisions
                if rand < BP %Binding probability
                    break;
                end
            end
        end
        SPTs{m} = particle_trace;
        all_particle_no_step(m,:) = i;
        all_particle_no_collision(m,:) = n;
    end
    
    for m= 1:n_sim
        SPT= SPTs{m};
        plot(SPT(:,1), SPT(:,2));  % plot trajectories
        plot(SPT(end,1), SPT(end,2), 'xk', 'MarkerSize', 15); %plot end points
        all_particle_final_point(m,:) = SPT(end,:);
        captured = captured + (SPT(end,2) >= -L);
    end
toc;