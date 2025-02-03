function dataMat = arrayIndex(dataMat, idx, dim)
% index from any dimension 'dim' in 'dataIn', using the index 'idx'.
% 'idx' can be an enumerated or logical index and should have not more
% entries as the target dimension 'dim'.
% Usage: dataOut = arrayIndex(dataIn, idx, dim)

cSize = size(dataMat); %this is the size of the source matrix
if length(cSize) < dim
    error('dim is larger as the number of available dimensions in dataIn. Use different dimension for indexing.')
end

%% check if trials is logical index.
if ~islogical(idx)
    temp = false(1, max(idx));
    temp(idx) = true;
    idx = temp;
end

%% check that idx matches size of target dimension
if length(idx) ~= cSize(dim)
    error('index does not match the size of the target dimension dim. Adjust index or change target dimension.')
end
    
%% do the indexing
dataMat = reshape(dataMat, prod(cSize(1:dim-1)), cSize(dim), []); %this reshapes the data matrix so that the indexing dimension is ensured to be dim 2
cSize(dim) = sum(idx); %this is the size of the target matrix
dataMat = reshape(dataMat(:,idx,:),cSize); %isolate target trials and reshape into size of the target matrix


