function [data, col] = degap( data, col)
%Degap database

    if size(data, 1) > length(col)
        %we need to transpose matrix
        data = data';
    end

    %Start degupping
    while true
        mat = isnan(data);
        [n , m] = size(mat);
        % Calculate fraction of missed for records
        mr = sum(mat, 2) / m;
        % Calculate fraction of missed for features
        mf = sum(mat) / n;
        if sum(mr > 0) == 1
            ind = mr > 0;
            data(ind, :) = [];
            col(ind) = [];
            break;
        end
        if sum(mf > 0) == 1
            ind = mf > 0;
            data(:, ind) = [];
            break;
        end

        %Search maximally gapped
        [mrm, imr] = max(mr);
        [mfm, imf] = max(mf);

        fprintf('Rec. %g  Attr. %g  Max rec. %g  Max attr. %g\n', sum(mr > 0), sum(mf > 0), mrm, mfm);

        if mrm > mfm
            data(imr, :) = [];
            col(imr) = [];
        else
            data(:, imf) = [];
        end
    end
end