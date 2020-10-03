%[str] = getMovieName(data) returns the identifier string ' date movieName' for each movie in data

% Francois Aguet 08/2013

function str = getMovieName(data)

nd = numel(data);
str = cell(1,nd);
for i = 1:nd
    if isempty(data(i).date)
        str{i} = [' ' getCellDir(data(i))];
    else
        str{i} = [' ' num2str(data(i).date) filesep getCellDir(data(i))];
    end
end
if nd==1
    str = str{1};
end