function numOut = getNumFromNumeral(stringIn)

if isequal(stringIn, 'I')
    numOut = 1;
elseif isequal(stringIn, 'II')
    numOut = 2;
elseif isequal(stringIn, 'III')
    numOut = 3;
elseif isequal(stringIn, 'IV')
    numOut = 4;
elseif isequal(stringIn, 'V')
    numOut = 5;
elseif isequal(stringIn, 'VI')
    numOut = 6;
elseif isequal(stringIn, 'VII')
    numOut = 7;
elseif isequal(stringIn, 'VIII')
    numOut = 8;
else 
    numOut = NaN;
end