PrintModel := proc (x, u, f, g, h, xnull, unull) local i, n, m, printf1, 
printh1; printf("The model  Xdot = F(x) +G(x)*U        \n"); if 4 < _npassed 
then printf("             y   = %a(%a)             is:   \n",h,x); end if; n :=
Dimension(f); printf1 := Vector(n); for i to n do printf1[i] := f[i]; end do; 
print(F^cat(`<`,DEAstep,`>`) = f); print(G^cat(`<`,DEAstep,`>`) = g); if 4 < 
_npassed then m := Dimension(h); printh1 := Vector(m); for i to m do printh1[i]
:= h[i]; end do; printf("%a(%a) = ",h,x), print(printh1); end if; if 5 < 
_npassed then printf("the working point xnull is %a \n",xnull); end if; if 6 <
_npassed then if nonaffine = true then printf("the working point unull is %a \n\
",unull); end if; end if; end proc;
