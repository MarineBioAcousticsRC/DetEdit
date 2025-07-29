function batches = split_into_batches(idxs)
    idxs = sort(idxs); % Ensure indices are sorted
    n = numel(idxs);
    batches = {};

    i = 1;
    while i <= n
        j = i + 1;
        if j > n
            % Only one index left
            batches{end+1} = [idxs(i), 1, idxs(i)];
            break;
        end

        % Determine step size between current and next index
        step = idxs(j) - idxs(i);
        while j < n && (idxs(j+1) - idxs(j)) == step
            j = j + 1;
        end

        % Add batch: start:step:end
        batches{end+1} = [idxs(i), step, idxs(j)];
        i = j + 1;
    end
end