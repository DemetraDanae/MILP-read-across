function D=fastdist(x,y)

% Calculates squared Euclidean distances between two sets of objetcs

D=((sum(y'.^2))'*ones(1,size(x,1)))+(ones(size(y,1),1)*(sum(x'.^2)))-2*(y*x');
