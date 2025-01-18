function frameRecon = svdFrameReconstruct_GN(U, V)
% function frameRecon = svdFrameReconstruct(U, V)
% U is Y x X x nSVD
% V is nSVD x nFrames

% reshape U to be nPix x nSVD
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));

if length(size(V)) < 3
    % this is the regular case for reconstructing a sequence of frames or
    % just a single frame

    % multiply and reshape back into Y x X
    frameRecon = reshape(Ur*V, size(U,1), size(U,2), size(V,2));
elseif length(size(V)) == 3
    % This unpacks higher Dimensional compilations of frames. 
    % e.g.: [n_frames, modality, further-conditions]
    
    % I have only tested this for len(size(V)) == 3, but it should also
    % work for anythink higher
    Vr = reshape(V, size(V, 1), prod(size(V, 2:length(size(V)))));
    
    % multiply and reshape back into Y x X
    frameRecon = reshape(Ur*Vr, cat(2, size(U,1:2), size(V, 2:length(size(V)))));
else
    error('Unsupported Nr. of V Dimensions');
end
