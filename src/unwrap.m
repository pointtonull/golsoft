function unwrapped=unwrap(wrapped)
    unwrapped = unwrap_wls(wrapped);

%    % debug block
%    tic; unwrapped_ls = unwrap_ls(wrapped); lstime = toc;
%    lstime
%    tic;unwrapped_wls = unwrap_wls(wrapped); wlstime = toc;
%    wlstime
%    fprintf('Time ls / timewls = %f\n', lstime / wlstime)
%    figure; imshow(mat2gray(wrapped));
%    figure; imshow(mat2gray(unwrapped_ls));
%    figure; imshow(mat2gray(unwrapped_wls));
%    diference = sum(sum((unwrapped_ls - unwrapped_wls) .^ 2));
%    diference
%    unwrapped = unwrapped_wls;

end
