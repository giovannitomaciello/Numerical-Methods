function [q,p] = cellForest(dTdq,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,updateGhost,savingStep)

    [ND, NP] = size(ptcls.x);

    q = zeros(ND,NP,nTime/savingStep);
    p = zeros(ND,NP,nTime/savingStep);
     
    %Ftmp ha dimensioni [ND,NP] e Ftmp(i,j) rappresenta la componente i sulla particella j  
    Ftmp = - dTdq(ptcls, grd_to_ptcl);

    t = 0;

    gamma = 1/(2 - 2^(1/3));

    for i = 2:nTime
        %momentum tmp
        ptmp = ptcls.p(:,1:NP) + gamma*dt/2*Ftmp(:,1:NP);
        
        %update position
        ptcls.x(:,1:NP) = ptcls.x(:,1:NP) + gamma * dt * dKdp(ptmp);

        %update boundary conditions
        ptcls.x(:,1:NP) = boundaryConditions(ptcls.x(:,1:NP));

        %update ghost particles
        ptcls.x = updateGhost(ptcls.x(:,1:NP), NP);

        %recalculate the grid
        grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

        %Force tmp
        Ftmp = - dTdq(ptcls, grd_to_ptcl);
       
        %momentum
        ptmp = ptmp + (1 - gamma) * dt/2 * Ftmp(:,1:NP);

        %update position
        ptcls.x(:,1:NP) = ptcls.x(:,1:NP) + (1 - 2*gamma) * dt * dKdp(ptmp);

        %update boundary conditions
        ptcls.x(:,1:NP) = boundaryConditions(ptcls.x(:,1:NP));

        %update ghost particles
        ptcls.x = updateGhost(ptcls.x(:,1:NP), NP);

        %recalculate the grid
        grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

        %Force
        Ftmp = - dTdq(ptcls, grd_to_ptcl);

        %momentum
        ptmp = ptmp + (1 - gamma) * dt/2 * Ftmp(:,1:NP);

        %update position
        ptcls.x(:,1:NP) = ptcls.x(:,1:NP) + gamma * dt * dKdp(ptmp);

        %update boundary conditions
        ptcls.x(:,1:NP) = boundaryConditions(ptcls.x(:,1:NP));

        %update ghost particles
        ptcls.x = updateGhost(ptcls.x(:,1:NP), NP);

        %recalculate the grid
        grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

        %Force
        Ftmp = - dTdq(ptcls, grd_to_ptcl);

        %momentum
        ptcls.p(:,1:NP) = ptmp + gamma * dt/2 * Ftmp(:,1:NP);

        % saving the results
        if mod(i,savingStep) == 0
            q(:,:,i/savingStep) = ptcls.x(:,1:NP);
            p(:,:,i/savingStep) = ptcls.p(:,1:NP);

            disp("-------------------------------------")
            % print the time
            fprintf('Time: %f\n', t);
            disp("-------------------------------------")
        end
        % save the time
        t = t + dt;
    end
end