function [A_new] = removeNaNs(A_old)
  A_new = A_old;
  nan_idx = isnan(A_new);
  height = size(A_new, 1);
  rm_idx = zeros(height, 1);

  for i = 1 : height
    if ~isempty( find( nan_idx(i, :), 1 ) )
      rm_idx(i) = 1;
    end
  end

  A_new( logical(rm_idx), : ) = [];
end
