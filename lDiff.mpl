Ldiff := proc (f, h, x, symp) local n, j, ld; n := Dimension(f); if type(h,'
table') then if not Dimension(h) = 1 then ERROR(
"2nd argument must be a scalar valued function"); end if; end if; ld := 0; for
j to n do ld := ld+diff(h,x[j])*f[j]; end do; if symp = 1 then simplify(normal(
ld)); end if; end proc;
