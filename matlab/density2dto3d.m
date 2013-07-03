function [densityGrid] = density2dto3d(img)

    imgSize = size(img);
    
    densityGrid = zeros(imgSize(1),imgSize(2),imgSize(1));

    for x=1:imgSize(1)
        for y=1:imgSize(2)
            
            densityGrid(x,y,1) = 1 - (img(x,y,1) + img(x,y,2) + img(x,y,3))/(3 * 255);
            
        end
    end

end