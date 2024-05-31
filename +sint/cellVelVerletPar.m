function [q,p] = cellVelVerletPar(dTdq,homogdTdq,rCut2,dKdp,dt,nTime,grd,ptcls,grd_to_ptcl,boundaryConditions,updateGhost,savingStep,hf,np)

    [ND, NP] = size(ptcls.x);

    %ptcls = structfun(@gpuArray, ptcls, 'UniformOutput', false);

    q = zeros(ND,NP,nTime/savingStep);
    p = zeros(ND,NP,nTime/savingStep);
     
    %Ftmp ha dimensioni [ND,NP] e Ftmp(i,j) rappresenta la componente i sulla particella j  
    Ftmp = - sint.forceCells(dTdq, ptcls, grd_to_ptcl, rCut2, grd.removeIndex);
    if hf
        Ftmp = Ftmp - homogdTdq(ptcls);
    end

    [nLy, nLx, nLz] = size(grd_to_ptcl);

    % prepare the indices for the parallel computation
    IndexX = cell(np(1),1);
    IndexY = cell(np(2),1);
    IndexZ = cell(np(3),1);

    for i = 1:np(1)
        IndexX{i} = (i-1)*nLx/np(1):i*nLx/np(1)+1;
        if i == 1
            IndexX{i}(1) = [];
        elseif i == np(1)
            IndexX{i}(end) = [];
        end
    end

    for i = 1:np(2)
        IndexY{i} = (i-1)*nLy/np(2):i*nLy/np(2)+1;
        if i == 1
            IndexY{i}(1) = [];
        elseif i == np(2)
            IndexY{i}(end) = [];
        end
    end

    for i = 1:np(3)
        IndexZ{i} = (i-1)*nLz/np(3):i*nLz/np(3)+1;
        if i == 1
            IndexZ{i}(1) = [];
        elseif i == np(3)
            IndexZ{i}(end) = [];
        end
    end

    [indexToSPMDX, indexToSPMDY, indexToSPMDZ] = ind2sub(np,1:prod(np));

    t = 0;

    for i = 2:nTime
        %momentum tmp
        ptmp = Ftmp(:,1:NP)*dt/2 + ptcls.p(:,1:NP);
       
        % position
        ptcls.x(:,1:NP) = dKdp(ptmp)*dt + ptcls.x(:,1:NP);

        %update boundary conditions
        ptcls.x(:,1:NP) = boundaryConditions(ptcls.x(:,1:NP));

        %update ghost particles
        ptcls = updateGhost(ptcls, NP);

        %recalculate the grid
        grd_to_ptcl = sint.init_ptcl_mesh(grd, ptcls);

        %Force tmp
        spmd (prod(np))
            grd_to_ptclPar = grd_to_ptcl(IndexX{indexToSPMDX(spmdIndex)},IndexY{indexToSPMDY(spmdIndex)},IndexZ{indexToSPMDZ(spmdIndex)});
            Ftmp_spmd = - sint.forceCells(dTdq, ptcls, grd_to_ptclPar, rCut2, grd.removeIndex);
            Ftmp = spmdPlus(Ftmp_spmd);
        end

        Ftmp = Ftmp{1};
        
        if hf
            Ftmp = Ftmp - homogdTdq(ptcls);
        end
       
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