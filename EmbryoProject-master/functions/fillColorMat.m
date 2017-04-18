
function colorMatOut = fillColorMat(colorMatIn, groupNumI, pnListI, i)

switch groupNumI
    case 1
        if pnListI == 0
            colorMatIn(i,:) = [1 .6 .6];
        elseif pnListI == 1
            colorMatIn(i,:) = [.8 .3 .3];
        else
            colorMatIn(i,:) = [.6 0 0];
        end
    case 2
        if pnListI == 0
            colorMatIn(i,:) = [.6 .6 1];
        elseif pnListI == 1
            colorMatIn(i,:) = [.3 .3 .8];
        else
            colorMatIn(i,:) = [0 0 .6];
        end
    case 3
        if pnListI == 0
            colorMatIn(i,:) = [.6 1 .6];
        elseif pnListI == 1
            colorMatIn(i,:) = [.3 .8 .3];
        else
            colorMatIn(i,:) = [0 .6 0];
        end
    case 4
        if pnListI == 0
            colorMatIn(i,:) = [1 1 .6];
        elseif pnListI == 1
            colorMatIn(i,:) = [.8 .8 .3];
        elseif pnListI == 2
            colorMatIn(i,:) = [.6 .6 0];
        else
            colorMatIn(i,:) = [.5 .5 .5];
        end
end

colorMatOut = colorMatIn;


