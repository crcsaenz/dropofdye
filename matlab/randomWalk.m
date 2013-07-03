function [densityGrid] = randomWalk (nx, ny, nz, numParticles, numSteps)

    densityGrid = zeros(nx, ny, nz);
    
    % perform a random walk for each of numParticles
    for i=1:numParticles
        
        % start at top and center (for now)
        currX = floor(nx/2);
        currY = ny;
        currZ = floor(nz/2);
        currDirection = 5; % 0 = right(+x), 1 = left(-x), 2 = forward(+z),
                            % 3 = backward(-z), 4 = up(+y), 5 = down(-y)
        
        % single step for given particle
        for step=0:numSteps-1
            
            p1 = rand(1);
            if (p1 < 0.5)
                % don't move, reset direction
                currDirection = 5;
            elseif (p1 < 0.8)
                % continue with current direction
                switch currDirection
                    case 0
                        if (currX < nx)
                            currX = currX + 1;
                        end
                    case 1
                        if (currX > 1)
                            currX = currX - 1;
                        end
                    case 2
                        if (currZ < nz)
                            currZ = currZ + 1;
                        end
                    case 3
                        if (currZ > 1)
                            currZ = currZ - 1;
                        end
                    case 4
                        if (currY < ny)
                            currY = currY + 1;
                        end
                    case 5
                        if (currY > 1)
                            currY = currY - 1;
                        end
                end
            else
                % pick random new direction
                p2 = rand(1);
                if (p2 < 0.3)
                    if (currY < ny)
                        currY = currY + 1;
                        currDirection = 4;
                    end
                elseif (p2 < 0.4)
                    if (currZ < nz)
                        currZ = currZ + 1;
                        currDirection = 2;
                    end
                elseif (p2 < 0.5)
                    if (currZ > 1)
                        currZ = currZ - 1;
                        currDirection = 3;
                    end
                elseif (p2 < 0.6)
                    if (currX < nx)
                        currX = currX + 1;
                        currDirection = 0;
                    end
                elseif (p2 < 0.7)
                    if (currX > 1)
                        currX = currX - 1;
                        currDirection = 1;
                    end
                else
                    if (currY > 1)
                        currY = currY - 1;
                        currDirection = 5;
                    end
                end %end random p2
                
            end %end random p1
            
            % increment densityGrid at current location
            densityGrid(currX,currY,currZ) = densityGrid(currX,currY,currZ) + 1;
            
        end %end numSteps loop
        
    end %end numParticles loop
    
    % scale densityGrid
    %densityGrid = densityGrid ./ (nx*ny*nz / numParticles);
    
    fileID = fopen('../density tests/dyeDensityTest27.pbrt','w');
    fprintf(fileID,'Volume "volumegrid" "integer nx" %d "integer ny" %d "integer nz" %d "point p0" [ 0.510000 -0.490000 0.010000 ] "point p1" [ 1.290000 1.490000 0.790000 ] "float density" [',nx,ny,nz);

    for z=1:nz
        for y=1:ny
            for x=1:nx
                fprintf(fileID,'%6.4f ',densityGrid(x,y,z));
            end
        end
    end
    
    fprintf(fileID,']');
    fclose(fileID);
    
end