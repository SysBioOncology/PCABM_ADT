function [cells,prop] = shuffleCells(cells,prop)

shf = randperm(length(cells));   %Prepare random shuffling
cells = cells(shf);              %Randomly shuffle cells

%for each cell property, shuffle accordingly
propnames=fieldnames(prop);

for n = 1:length(propnames)
    prop.(propnames{n})=prop.(propnames{n})(shf);
end

end
