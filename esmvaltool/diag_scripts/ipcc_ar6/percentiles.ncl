function weighted_mean(x, w)
  local x, w, y
begin
  y = dim_sum_n(x*conform(x, w, 0), 0)
  y = y / dim_sum_n(where(ismissing(x), w@_FillValue, conform(x, w, 0)), 0)
  return(y)
end

function weighted_percentile(x, w, t)
local x, w, xx, ww, wsum, j, k, dims, idx
begin
  xx = x
  dims = dimsizes(xx)
  if (dimsizes(dims).eq.1) then
    idx = dim_pqsort(xx, 2)
    wsum = sum(where(ismissing(x), w@_FillValue, w))
    ww = 0.
    do j = 0, dims(0)-1
      if (.not.ismissing(xx(j))) then
        ww = ww + w(idx(j))
        if (ww/wsum.eq.t) then
          y = xx(j)
          break
        elseif (ww/wsum.gt.t) then
          if (j.eq.0) then
            y = xx(j)
          else
            y = (xx(j) + xx(j-1))/2.
          end if
          break
        end if
      end if
    end do
  else
    y = new(dims(1), typeof(x))
    idx = dim_pqsort_n(xx, 2, 0)
    do k = 0, dims(1)-1
      wsum = sum(where(ismissing(x(:, k)), w@_FillValue, w))
      ww = 0.
      do j = 0, dims(0)-1
        if (.not.ismissing(xx(j, k))) then
          ww = ww + w(idx(j, k))
          if (ww/wsum.eq.t) then
            y(k) = xx(j, k)
            break
          elseif (ww/wsum.gt.t) then
            if (j .gt. 0) then
              y(k) = (xx(j, k) + xx(j-1, k))/2.
            else
              y(k) = xx(j, k)
            end if
            break
          end if
        end if
      end do
    end do
  end if
  return(y)
end
