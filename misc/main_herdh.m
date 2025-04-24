% Hybrid Encrypted Reversible Data Hiding (HERDH) Implementation in MATLAB

% Main function to demonstrate the HERDH scheme
function main_herdh()
    % Clear workspace and command window
    clear all;
    close all;
    clc;
    
    % Load test images
    pentagon = imread('/MATLAB Drive/pentagon.jpeg');
    barbara = imread('/MATLAB Drive/barbara_gray.bmp');
    boat = imread('/MATLAB Drive/boat.jpeg');
    lena = imread('/MATLAB Drive/lena.jpg');
    
    % Convert to grayscale if needed and resize to standard size
    pentagon = prepare_image(pentagon);
    barbara = prepare_image(barbara);
    boat = prepare_image(boat);
    lena = prepare_image(lena);
    
    % Create confidential data (grayscale)
    confidential_gray = uint8(zeros(256, 256));
    for i = 1:256
        for j = 1:256
            confidential_gray(i, j) = mod((i * 3) + (j * 2), 256);
        end
    end
    
    % Create confidential data (color)
    confidential_color = uint8(zeros(256, 256, 3));
    for i = 1:256
        for j = 1:256
            confidential_color(i, j, 1) = mod(i * 2, 256);  % Red
            confidential_color(i, j, 2) = mod(j * 2, 256);  % Green
            confidential_color(i, j, 3) = mod((i + j) * 2, 256);  % Blue
        end
    end
    
    % Process each cover image
    process_image(pentagon, confidential_gray, confidential_color, 'pentagon');
    process_image(barbara, confidential_gray, confidential_color, 'barbara');
    process_image(boat, confidential_gray, confidential_color, 'boat');
    process_image(lena, confidential_gray, confidential_color, 'lena');
    
    % Create comparison tables and graphs
    create_comparison_tables();
end

% Function to prepare an image (convert to grayscale and resize)
function img = prepare_image(img)
    % Convert to grayscale if it's a color image
    if size(img, 3) > 1
        img = rgb2gray(img);
    end
    
    % Resize to standard size (512x512)
    img = imresize(img, [512, 512]);
end

% Function to process an image with HERDH
% Update the process_image function to fix the overlapping titles
function process_image(cover_image, confidential_gray, confidential_color, name)
    % Process with grayscale confidential data
    [embedded_gray, svd_components_gray] = hide_data(cover_image, confidential_gray);
    [recovered_ci_gray, recovered_cdi_gray] = extract_data(embedded_gray, svd_components_gray, size(confidential_gray), false);
    
    % Calculate metrics for grayscale
    mse_gray = calculate_mse(cover_image, embedded_gray);
    psnr_embedded_gray = calculate_psnr(cover_image, embedded_gray);
    psnr_restored_gray = calculate_psnr(cover_image, recovered_ci_gray);
    entropy_original = calculate_entropy(cover_image);
    entropy_stego_gray = calculate_entropy(embedded_gray);
    
    % Process with color confidential data
    [embedded_color, svd_components_color] = hide_data(cover_image, confidential_color);
    [recovered_ci_color, recovered_cdi_color] = extract_data(embedded_color, svd_components_color, [size(confidential_color, 1), size(confidential_color, 2)], true);
    
    % Calculate metrics for color
    mse_color = calculate_mse(cover_image, embedded_color);
    psnr_embedded_color = calculate_psnr(cover_image, embedded_color);
    psnr_restored_color = calculate_psnr(cover_image, recovered_ci_color);
    entropy_stego_color = calculate_entropy(embedded_color);
    
    % Store results in a file
    results = struct('name', name, ...
                    'mse_gray', mse_gray, ...
                    'mse_color', mse_color, ...
                    'psnr_embedded_gray', psnr_embedded_gray, ...
                    'psnr_embedded_color', psnr_embedded_color, ...
                    'psnr_restored_gray', psnr_restored_gray, ...
                    'psnr_restored_color', psnr_restored_color, ...
                    'entropy_original', entropy_original, ...
                    'entropy_stego_gray', entropy_stego_gray, ...
                    'entropy_stego_color', entropy_stego_color);
    
    save([name '_results.mat'], 'results');
    
    % Create a figure with better spacing
    figure('Name', [name ' - Results'], 'Position', [100, 100, 1200, 800]);
    
    % Add a main title for grayscale results
    subplot(2, 1, 1);
    text(0.5, 0.95, 'Results with Grayscale Confidential Data', 'FontSize', 14, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off;
    
    % Original and embedded images (gray confidential data)
    subplot(2, 4, 1);
    imshow(cover_image);
    title('Original Image', 'FontSize', 10);
    
    subplot(2, 4, 2);
    imshow(embedded_gray);
    title({'Embedded Image'; ['PSNR: ' num2str(psnr_embedded_gray, '%.2f') 'dB']}, 'FontSize', 10);
    
    subplot(2, 4, 3);
    imshow(recovered_ci_gray);
    title({'Recovered Image'; ['PSNR: ' num2str(psnr_restored_gray, '%.2f') 'dB']}, 'FontSize', 10);
    
    subplot(2, 4, 4);
    imshow(confidential_gray);
    title('Confidential Data', 'FontSize', 10);
    
    % Add a main title for color results
    subplot(2, 1, 2);
    text(0.5, 0.95, 'Results with Color Confidential Data', 'FontSize', 14, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Units', 'normalized');
    axis off;
    
    % Original and embedded images (color confidential data)
    subplot(2, 4, 5);
    imshow(cover_image);
    title('Original Image', 'FontSize', 10);
    
    subplot(2, 4, 6);
    imshow(embedded_color);
    title({'Embedded Image'; ['PSNR: ' num2str(psnr_embedded_color, '%.2f') 'dB']}, 'FontSize', 10);
    
    subplot(2, 4, 7);
    imshow(recovered_ci_color);
    title({'Recovered Image'; ['PSNR: ' num2str(psnr_restored_color, '%.2f') 'dB']}, 'FontSize', 10);
    
    subplot(2, 4, 8);
    imshow(confidential_color);
    title('Confidential Data', 'FontSize', 10);
    
    % Adjust spacing between subplots
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.1, 0.1, 0.8, 0.8]);
    
    % Save figure
    saveas(gcf, [name '_results.fig']);
    saveas(gcf, [name '_results.png']);
    
    % Create histogram
    figure('Name', [name ' - Histogram']);
    imhist(cover_image);
    title(['Histogram of ' name]);
    saveas(gcf, [name '_histogram.fig']);
    saveas(gcf, [name '_histogram.png']);
end

% Function to encrypt an image
function encrypted = encrypt_image(image)
    % Simple XOR encryption with a key
    key = 123;
    encrypted = image;
    
    if ndims(image) == 2  % Grayscale
        [rows, cols] = size(encrypted);
        for i = 1:rows
            for j = 1:cols
                encrypted(i, j) = bitxor(encrypted(i, j), key);
            end
        end
    else  % Color
        [rows, cols, channels] = size(encrypted);
        for i = 1:rows
            for j = 1:cols
                for k = 1:channels
                    encrypted(i, j, k) = bitxor(encrypted(i, j, k), key);
                end
            end
        end
    end
end

% Function to decrypt an image
function decrypted = decrypt_image(encrypted_image)
    % XOR is its own inverse with the same key
    decrypted = encrypt_image(encrypted_image);
end

% Function to hide data in a cover image
function [embedded_image, svd_components] = hide_data(cover_image, confidential_data)
    % Convert images to double for processing
    cover_image_double = double(cover_image);
    
    % Step 1: Apply DWT over Cover Image
    [cA, cH, cV, cD] = dwt2(cover_image_double, 'haar');
    
    % Step 2: Encrypt Confidential Data Image
    encrypted_cdi = encrypt_image(confidential_data);
    
    % Convert to double after encryption
    if ndims(encrypted_cdi) == 3  % Color image
        % Convert to grayscale for embedding
        encrypted_cdi_gray = double(rgb2gray(encrypted_cdi));
    else  % Grayscale
        encrypted_cdi_gray = double(encrypted_cdi);
    end
    
    % Step 3: Apply SVD over HH band of Cover Image
    [U_hh, S_hh, V_hh] = svd(cD);
    
    % Step 4: Apply SVD over encrypted confidential data
    % Resize encrypted_cdi to match HH dimensions if needed
    if size(encrypted_cdi_gray) ~= size(cD)
        encrypted_cdi_gray = imresize(encrypted_cdi_gray, size(cD));
    end
    
    [U_cdi, S_cdi, V_cdi] = svd(encrypted_cdi_gray);
    
    % Get diagonal values from S matrices
    s_hh_diag = diag(S_hh);
    s_cdi_diag = diag(S_cdi);
    
    % Adjust S_cdi length to match S_hh if needed
    min_length = min(length(s_hh_diag), length(s_cdi_diag));
    s_hh_adj = s_hh_diag(1:min_length);
    s_cdi_adj = s_cdi_diag(1:min_length);
    
    % Step 5: Add SVD matrix of CDI and ECDI
    % We're adding the singular values only
    alpha = 0.01;  % Scaling factor to control embedding strength
    s_combined = s_hh_adj + alpha * s_cdi_adj;
    
    % Create new S matrix
    S_combined = zeros(size(S_hh));
    S_combined(1:min_length, 1:min_length) = diag(s_combined);
    
    % Step 6: Apply inverse SVD over combined SVD matrix
    HH_modified = U_hh * S_combined * V_hh';
    
    % Step 7: Apply inverse DWT to generate embedded Cover Image
    embedded_image_double = idwt2(cA, cH, cV, HH_modified, 'haar', size(cover_image_double));
    
    % Convert back to uint8
    embedded_image = uint8(embedded_image_double);
    
    % Return embedded image and SVD components needed for extraction
    svd_components = struct('U_hh', U_hh, 'S_hh', S_hh, 'V_hh', V_hh, ...
                           'min_length', min_length, 'alpha', alpha, ...
                           'U_cdi', U_cdi, 'V_cdi', V_cdi);
end

% Function to extract data from an embedded image
function [recovered_ci, recovered_cdi] = extract_data(embedded_image, svd_components, confidential_data_shape, is_color)
    % Unpack SVD components
    U_hh_orig = svd_components.U_hh;
    S_hh_orig = svd_components.S_hh;
    V_hh_orig = svd_components.V_hh;
    min_length = svd_components.min_length;
    alpha = svd_components.alpha;
    U_cdi = svd_components.U_cdi;
    V_cdi = svd_components.V_cdi;
    
    % Convert to double for processing
    embedded_image_double = double(embedded_image);
    
    % Step 1: Apply DWT over embedded image
    [cA, cH, cV, cD] = dwt2(embedded_image_double, 'haar');
    
    % Step 2: Apply SVD over HH band of embedded image
    [U_hh_emb, S_hh_emb, V_hh_emb] = svd(cD);
    
    % Get diagonal values
    s_hh_orig_diag = diag(S_hh_orig);
    s_hh_emb_diag = diag(S_hh_emb);
    
    % Step 3: Subtract SVD matrix of CI from embedded image
    s_hh_orig_adj = s_hh_orig_diag(1:min_length);
    s_hh_emb_adj = s_hh_emb_diag(1:min_length);
    
    s_cdi = (s_hh_emb_adj - s_hh_orig_adj) / alpha;
    
    % Create S matrix for original HH band
    S_hh_orig_full = zeros(size(S_hh_orig));
    S_hh_orig_full(1:length(s_hh_orig_diag), 1:length(s_hh_orig_diag)) = diag(s_hh_orig_diag);
    
    % Step 4: Reconstruct original HH band
    HH_orig = U_hh_orig * S_hh_orig_full * V_hh_orig';
    
    % Step 5: Apply inverse DWT to recover original cover image
    recovered_ci_double = idwt2(cA, cH, cV, HH_orig, 'haar', size(embedded_image_double));
    recovered_ci = uint8(recovered_ci_double);
    
    % Step 6: Reconstruct encrypted confidential data
    % Create S matrix for confidential data
    S_cdi = zeros(size(cD));
    S_cdi(1:min_length, 1:min_length) = diag(s_cdi);
    
    encrypted_cdi_double = U_cdi * S_cdi * V_cdi';
    
    % Resize to match the original confidential data shape
    encrypted_cdi_double = imresize(encrypted_cdi_double, confidential_data_shape);
    encrypted_cdi = uint8(encrypted_cdi_double);
    
    if is_color
        % For color images, we need to reconstruct the color channels
        % This is a simplified approach - in practice, you might need a more sophisticated method
        encrypted_cdi = repmat(encrypted_cdi, [1, 1, 3]);
    end
    
    % Step 7: Decrypt to get confidential data
    recovered_cdi = decrypt_image(encrypted_cdi);
end

% Function to calculate MSE between two images
function mse = calculate_mse(original, modified)
    % Ensure both images are grayscale
    if size(original, 3) > 1
        original = rgb2gray(original);
    end
    if size(modified, 3) > 1
        modified = rgb2gray(modified);
    end
    
    % Ensure both images are the same size
    if ~isequal(size(original), size(modified))
        modified = imresize(modified, size(original));
    end
    
    original_double = double(original);
    modified_double = double(modified);
    mse = mean((original_double(:) - modified_double(:)).^2);
end

% Function to calculate PSNR between two images
function psnr_val = calculate_psnr(original, modified)
    % Ensure both images are grayscale
    if size(original, 3) > 1
        original = rgb2gray(original);
    end
    if size(modified, 3) > 1
        modified = rgb2gray(modified);
    end
    
    % Ensure both images are the same size
    if ~isequal(size(original), size(modified))
        modified = imresize(modified, size(original));
    end
    
    mse = calculate_mse(original, modified);
    if mse == 0
        psnr_val = 100;  % Arbitrary high value for identical images
    else
        max_pixel = 255;
        psnr_val = 10 * log10((max_pixel^2) / mse);
    end
end

% Function to calculate entropy of an image
function ent = calculate_entropy(image)
    if ndims(image) == 3  % Color image
        % Convert to grayscale for entropy calculation
        gray = rgb2gray(image);
    else
        gray = image;
    end
    
    [counts, ~] = imhist(gray);
    counts = counts / sum(counts);
    counts = counts(counts > 0);
    ent = -sum(counts .* log2(counts));
end

% Function to create comparison tables and graphs
function create_comparison_tables()
    % Load results
    pentagon_results = load('pentagon_results.mat');
    barbara_results = load('barbara_results.mat');
    boat_results = load('boat_results.mat');
    lena_results = load('lena_results.mat');
    
    % Existing method results (from the paper)
    mse_existing = [0.89, 1.38, 1.21, 1.07];
    psnr_embedded_existing = [10.08, 8.68, 9.43, 9.93];
    psnr_restored_existing = [49.65, 49.76, 49.86, 49.99];
    
    % Proposed method results
    mse_proposed = [pentagon_results.results.mse_gray, barbara_results.results.mse_gray, ...
                   boat_results.results.mse_gray, lena_results.results.mse_gray];
    
    psnr_embedded_proposed = [pentagon_results.results.psnr_embedded_gray, ...
                             barbara_results.results.psnr_embedded_gray, ...
                             boat_results.results.psnr_embedded_gray, ...
                             lena_results.results.psnr_embedded_gray];
    
    psnr_restored_proposed = [pentagon_results.results.psnr_restored_gray, ...
                             barbara_results.results.psnr_restored_gray, ...
                             boat_results.results.psnr_restored_gray, ...
                             lena_results.results.psnr_restored_gray];
    
    % Create MSE comparison table and graph
    figure('Name', 'MSE Comparison');
    bar([mse_existing; mse_proposed]');
    title('MSE Comparison');
    xlabel('Images');
    ylabel('MSE');
    set(gca, 'XTickLabel', {'Pentagon', 'Barbara', 'Boat', 'Lena'});
    legend('Existing Method', 'Proposed Method');
    saveas(gcf, 'mse_comparison.fig');
    saveas(gcf, 'mse_comparison.png');
    
    % Create PSNR embedded comparison table and graph
    figure('Name', 'PSNR of Encrypted Embedded Images');
    bar([psnr_embedded_existing; psnr_embedded_proposed]');
    title('PSNR of Encrypted Embedded Images');
    xlabel('Images');
    ylabel('PSNR (dB)');
    set(gca, 'XTickLabel', {'Pentagon', 'Barbara', 'Boat', 'Lena'});
    legend('Existing Method', 'Proposed Method');
    saveas(gcf, 'psnr_embedded_comparison.fig');
    saveas(gcf, 'psnr_embedded_comparison.png');
    
    % Create PSNR restored comparison table and graph
    figure('Name', 'PSNR of Restored Images');
    bar([psnr_restored_existing; psnr_restored_proposed]');
    title('PSNR of Restored Images');
    xlabel('Images');
    ylabel('PSNR (dB)');
    set(gca, 'XTickLabel', {'Pentagon', 'Barbara', 'Boat', 'Lena'});
    legend('Existing Method', 'Proposed Method');
    saveas(gcf, 'psnr_restored_comparison.fig');
    saveas(gcf, 'psnr_restored_comparison.png');
    
    % Create entropy comparison tables and graphs
    entropy_original = [pentagon_results.results.entropy_original, ...
                       barbara_results.results.entropy_original, ...
                       boat_results.results.entropy_original, ...
                       lena_results.results.entropy_original];
    
    entropy_stego_gray = [pentagon_results.results.entropy_stego_gray, ...
                         barbara_results.results.entropy_stego_gray, ...
                         boat_results.results.entropy_stego_gray, ...
                         lena_results.results.entropy_stego_gray];
    
    entropy_stego_color = [pentagon_results.results.entropy_stego_color, ...
                          barbara_results.results.entropy_stego_color, ...
                          boat_results.results.entropy_stego_color, ...
                          lena_results.results.entropy_stego_color];
    
    % Entropy comparison (Gray)
    figure('Name', 'Entropy Comparison (Gray)');
    bar([entropy_original; entropy_stego_gray]');
    title('Entropy Comparison (Gray Confidential Data)');
    xlabel('Images');
    ylabel('Entropy');
    set(gca, 'XTickLabel', {'Pentagon', 'Barbara', 'Boat', 'Lena'});
    legend('Original Image', 'Stego Image');
    saveas(gcf, 'entropy_gray_comparison.fig');
    saveas(gcf, 'entropy_gray_comparison.png');
    
    % Entropy comparison (Color)
    figure('Name', 'Entropy Comparison (Color)');
    bar([entropy_original; entropy_stego_color]');
    title('Entropy Comparison (Color Confidential Data)');
    xlabel('Images');
    ylabel('Entropy');
    set(gca, 'XTickLabel', {'Pentagon', 'Barbara', 'Boat', 'Lena'});
    legend('Original Image', 'Stego Image');
    saveas(gcf, 'entropy_color_comparison.fig');
    saveas(gcf, 'entropy_color_comparison.png');
    
    % Create comprehensive results table
    comprehensive_results = table(...
        {'Pentagon'; 'Barbara'; 'Boat'; 'Lena'}, ...
        mse_existing', mse_proposed', ...
        psnr_embedded_existing', psnr_embedded_proposed', ...
        psnr_restored_existing', psnr_restored_proposed', ...
        entropy_original', entropy_stego_gray', entropy_stego_color', ...
        'VariableNames', {'Image', 'MSE_Existing', 'MSE_Proposed', ...
                         'PSNR_Embedded_Existing', 'PSNR_Embedded_Proposed', ...
                         'PSNR_Restored_Existing', 'PSNR_Restored_Proposed', ...
                         'Entropy_Original', 'Entropy_Stego_Gray', 'Entropy_Stego_Color'});
    
    % Display comprehensive results
    disp(comprehensive_results);
    
    % Save comprehensive results
    writetable(comprehensive_results, 'comprehensive_results.csv');
    save('comprehensive_results.mat', 'comprehensive_results');
end