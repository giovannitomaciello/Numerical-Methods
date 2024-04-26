function [q,p] = cellVelVerlet(dTdq,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,updateGhost,savingStep)

    [ND, NP] = size(ptcls.x);

    q = zeros(ND,NP,nTime/savingStep);
    p = zeros(ND,NP,nTime/savingStep);
     
    %Ftmp ha dimensioni [ND,NP] e Ftmp(i,j) rappresenta la componente i sulla particella j  
    Ftmp = - dTdq(ptcls, grd_to_ptcl);

    t = 0;

    for i = 2:nTime
        %momentum tmp
        ptmp = Ftmp(:,1:NP)*dt/2 + ptcls.p(:,1:NP);
       
        % position
        ptcls.x(:,1:NP) = dKdp(ptmp)*dt + ptcls.x(:,1:NP);

        %update boundary conditions
        ptcls.x(:,1:NP) = boundaryConditions(ptcls.x(:,1:NP));

        %update ghost particles
        ptcls.x = updateGhost(ptcls.x(:,1:NP), NP);

        %recalculate the grid
        grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

        %Force tmp
        Ftmp = - dTdq(ptcls, grd_to_ptcl);
       
        %momentum
        ptcls.p(:,1:NP) = Ftmp(:,1:NP)*dt/2 + ptmp(:,1:NP);

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