% Read the picture of the tiger, and convert to black and white.
img = imread('image_compression/tiger.jpg') ; 

epsilon = 0.05;
norm_val = 2;

% compressed_img = compress_rgb(img, epsilon, norm_val);
compressed_img = compress_gray(img, epsilon, norm_val);

% epsilon is the error
function compressed_img = compress_rgb(img, epsilon, norm_val)
    img = im2double(img);
    [nx,ny,n] = size(img) ;
    compressed_img = zeros(nx,ny,n);
    
    R = img(:,:,1) ; 
    G = img(:,:,2) ; 
    B = img(:,:,3) ; 
    
    [ur,sr,vr]=svd(R) ;
    [ug,sg,vg]=svd(G) ;
    [ub,sb,vb]=svd(B) ;

    % Plot the magnitude of the singular values (log scale)
    sigmas_r = diag(sr);
    sigmas_g = diag(sg);
    sigmas_b = diag(sb);
    figure; plot(log10(sigmas_r)); title('Singular Values (red) (Log10 Scale)');
    figure; plot(log10(sigmas_g)); title('Singular Values (green) (Log10 Scale)');
    figure; plot(log10(sigmas_b)); title('Singular Values (blue) (Log10 Scale)');
    
    figure; plot(cumsum(sigmas_r) / sum(sigmas_r)); title('Cumulative Percent of Total Sigmas (red)');
    figure; plot(cumsum(sigmas_g) / sum(sigmas_g)); title('Cumulative Percent of Total Sigmas (green)');
    figure; plot(cumsum(sigmas_b) / sum(sigmas_b)); title('Cumulative Percent of Total Sigmas (blue)');
    
    % Show full-rank img
    figure; imshow(img), title('Full-Rank img');

    % Compute low-rank approximations of the img, and show them
    ns = length(sigmas_r);
    precision = zeros(ns);
    norm_r = norm(R, norm_val);
    norm_g = norm(G, norm_val);
    norm_b = norm(B, norm_val);
    for i = 1:ns
        % RED
        % Keep largest singular values, and nullify others.
        approx_sigmas_r = sigmas_r; 
        approx_sigmas_r(i+1:end) = 0;

        % Form the singular value matrix, padded as necessary
        approx_sr = sr; 
        approx_sr(1:ns, 1:ns) = diag(approx_sigmas_r);

        % Compute low-rank approximation by multiplying out component matrices.
        approx_img_r = ur * approx_sr * vr';
        
        % GREEN
        % Keep largest singular values, and nullify others.
        approx_sigmas_g = sigmas_g; 
        approx_sigmas_g(i+1:end) = 0;

        % Form the singular value matrix, padded as necessary
        approx_sg = sg; 
        approx_sg(1:ns, 1:ns) = diag(approx_sigmas_g);

        % Compute low-rank approximation by multiplying out component matrices.
        approx_img_g = ug * approx_sg * vg';
        
        % BLUE
        % Keep largest singular values, and nullify others.
        approx_sigmas_b = sigmas_b; 
        approx_sigmas_b(i+1:end) = 0;

        % Form the singular value matrix, padded as necessary
        approx_sb = sb; 
        approx_sb(1:ns, 1:ns) = diag(approx_sigmas_b);

        % Compute low-rank approximation by multiplying out component matrices.
        approx_img_b = ub * approx_sb * vb';

        
        compressed_img(:, :, 1) = approx_img_r;
        compressed_img(:, :, 2) = approx_img_g;
        compressed_img(:, :, 3) = approx_img_b;
        
        rel_err_r = norm(R - approx_img_r, norm_val) / norm_r;
        rel_err_g = norm(G - approx_img_g, norm_val) / norm_g;
        rel_err_b = norm(B - approx_img_b, norm_val) / norm_b;
        rel_err = (rel_err_r + rel_err_g + rel_err_b) / 3;
        precision(i) = 100 - rel_err * 100;
        
        if rel_err < epsilon
            % Plot approximation
            figure; imshow(compressed_img); title(sprintf('Rank %d img', i));
            break;
        end
    end
    
    figure; plot(precision); title('Relative Error Percentage');
end


function compressed_img = compress_gray(img, epsilon, norm_val)
    img = rgb2gray(img);
    img = im2double(img);
    [nx,ny] = size(img) ;
    compressed_img = zeros(nx,ny);

    % Compute SVD of this img
    [U, S, V] = svd(img);

    % Plot the magnitude of the singular values (log scale)
    sigmas = diag(S);
%     figure; plot(log10(sigmas)); title('Singular Values (Log10 Scale)');
%     
%     figure; plot(cumsum(sigmas) / sum(sigmas)); title('Cumulative Percent of Total Sigmas');
%     
%     % Show full-rank img
%     figure; imshow(img), title('Full-Rank img');

    % Compute low-rank approximations of the img, and show them
    ns = length(sigmas);    
    
    precision = zeros(ns);
    img_norm = norm(img, norm_val);
    for i = 1:ns
        % Keep largest singular values, and nullify others.
        approx_sigmas = sigmas; 
        approx_sigmas(i+1:end) = 0;

        % Form the singular value matrix, padded as necessary
        approx_S = S; 
        approx_S(1:ns, 1:ns) = diag(approx_sigmas);

        % Compute low-rank approximation by multiplying out component matrices.
        compressed_img = U * approx_S * V';
        
        rel_err = norm(img - compressed_img, norm_val) / img_norm;
        precision(i) = 100 - rel_err * 100;

        if rel_err < epsilon
            % Plot approximation
            figure; imshow(compressed_img); title(sprintf('Rank %d img', i));
            break;
        end
    end
    figure; plot(precision); title('Relative Error Percentage');
end
