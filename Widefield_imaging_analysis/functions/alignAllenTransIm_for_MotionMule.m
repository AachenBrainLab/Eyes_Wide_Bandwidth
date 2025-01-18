function im = alignAllenTransIm_for_MotionMule(im, transParams, targSize)
% im = alignAllenTransIm(im, transParams, targSize)
% 
% Take a brain image and rotate, scale, and translate it according to the
% parameters in transParams (produced by alignBrainToAllen GUI). im may be
% an image stack. Note that due to the behavior of imrotate, NaNs may not
% propagate identically if an image stack is used vs. calling this function
% on individual images.
%
% targSize defines the dimensions of the output image after alignment. This 
% is needed to match the size of the allen brain map, which is usually
% 540x586 pixels (default if targSize is not assigned).
%
% Pixels where the value is not defined (due to rotation) are set to NaN.
% This behavior is different from imrotate.

if ~exist('targSize', 'var') || isempty(targSize)
    targSize = [540, 586]; %standard size of allen brain map
end

offset = 5E1;
[h, w, d] = size(im);

%% This was removed on github ~2024-02, but it broke all the alignments!
% add NaNs if image is smaller as targSize
if h < targSize(1)
    im = [im; NaN(targSize(1) - h, w, d)];
    [h, w, d] = size(im);
end

if w < targSize(2)
    im = [im, NaN(h, targSize(2) - w, d)];
end


%%
% % add NaNs if image is smaller as targSize
% if h < targSize(1)
% %     im = [NaN(floor(targSize(1) - h)/2, w, d); im; NaN(floor(targSize(1) - h)/2, w, d)];
%     im = [im; NaN(floor(targSize(1) - h), w, d)];
%     [h, w, d] = size(im);
% end
% 
% if w < targSize(2)
%     im = [im, NaN(h, floor(targSize(1) - w), d)];
% %     im = [NaN(h, floor(targSize(1) - w)/2, d), im, NaN(h, floor(targSize(1) - w)/2, d)];
% end
dSize = size(im);

% Set pixels off from zero
[theMin, minIdx] = min(im(:));
im = im - theMin + offset;

% Rotate
im = imrotate(im, transParams.angleD, 'bilinear');

% Scale
if transParams.scaleConst ~= 1
  im = imresize(im, transParams.scaleConst);
end

% Set NaNs to 0, because imtranslate can't handle them
nans = isnan(im);
if any(nans(:))
  im(nans) = 0;
end

% Translate
im = imtranslate(im, transParams.tC');

% Detect missing pixels due to rotation, set to NaN
im(im < offset) = NaN;

% remove some pixels at the edges to avoid rotation artifact
im = removeOutline(im,10,0);

% Restore offset
im = im + theMin - offset;
im(minIdx) = theMin;

% Check if image needs to be trimed because rotate expands it or expanded
% because picture is smaller than expected.
newSize = size(im);
trimH = floor((newSize(1) - dSize(1)) / 2);
trimW = floor((newSize(2) - dSize(2)) / 2);

if trimH > 0
    im = im(trimH + (1:dSize(1)), :, :);
elseif trimH < 0
    im = cat(1, nan(-trimH, size(im,2), size(im,3)), im);
end

if trimW > 0
    im = im(:, trimW + (1:dSize(2)), :);
elseif trimW < 0
    im = cat(2, nan(size(im,1), -trimW , size(im,3)), im);
end

%crop image if larger as targSize
if size(im,1) > targSize(1)
    im = im(1:targSize(1), :, :);
end
if size(im,2) > targSize(2)
    im = im(:, 1:targSize(2), :);
end

% add nans if image got smaller 
if size(im,1) < targSize(1) || size(im,2) < targSize(2)
    nIm = nan(targSize(1), targSize(2), size(im,3), class(im));
    nIm(1 : size(im,1), 1 : size(im,2), :) = im;
    im = nIm;
end