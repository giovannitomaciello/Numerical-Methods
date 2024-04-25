function [q,p] = cellVelVerlet(dTdq,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,savingStep)

    [ND, NP] = size(ptcls.x);

    q = zeros(ND,NP,nTime/savingStep);
    p = zeros(ND,NP,nTime/savingStep);
     
    %Ftmp ha dimensioni [ND,NP] e Ftmp(i,j) rappresenta la componente i sulla particella j  
    Ftmp = - dTdq(ptcls, grd_to_ptcl);

    t = 0;

    for i = 2:nTime
        %momentum tmp
        ptmp = Ftmp*dt/2 + ptcls.p;
       
        % position
        ptcls.x = dKdp(ptmp)*dt + ptcls.x;
        
        %update boundary conditions
        ptcls = boundaryConditions(ptcls);

        %recalculate the grid
        grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

        %Force tmp
        Ftmp = - dTdq(ptcls, grd_to_ptcl);
       
        %momentum
        ptcls.p = Ftmp*dt/2 + ptmp;

        % saving the results
        if mod(i,savingStep) == 0
            q(:,:,i/savingStep) = ptcls.x;
            p(:,:,i/savingStep) = ptcls.p;

            disp("-------------------------------------")
            % print the time
            fprintf('Time: %f\n', t);
            disp("-------------------------------------")
        end
        % save the time
        t = t + dt;
    end
end