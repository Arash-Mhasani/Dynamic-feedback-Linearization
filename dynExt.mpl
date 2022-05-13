DynExt := proc (f, g, h, x, xnull, fn, gn, hn, xn, xnulln, Adegn, rdegn) local
n, m, y1, p, bdeg, rdeg, rtot, i, q1, nq, mq, j, kk, Acal, bcal, delta, alpha,
beta, fhelp, g1, g2, Adeghelp, xnhelp, g11, g12, g21, g22, q, qtest, pivot, 
beta1, beta2, fnew, gnew, k, xhelp, locxnulln, Adeg, dum, locfn, locgn, nrem, 
var1seq, reldegree, u, unull, xnew, K, S, K1, S1, io1, jo1, maplet1, maplet2, 
result, IntMaplet2, tmpU, SearchText, tmpResult, jj; global io, jo, DEAstep, 
numofNonZe; tmpResult := 0; maplet1 := Maplets:-Elements:-Maplet(Maplets:-
Elements:-InputDialog["ID1"]("Enter the expression for k(x), d=default (k__d=-L\
__(f^(itr-1))^(r^(itr-1)) h__(i0)):",(':-onapprove') = Maplets:-Elements:-
Shutdown(["ID1"]),(':-oncancel') = Maplets:-Elements:-Shutdown())); maplet2 :=
Maplets:-Elements:-Maplet(Maplets:-Elements:-InputDialog["ID1"](
"Enter the expression for s(x), d=default (`s__d`=1):",(':-onapprove') = 
Maplets:-Elements:-Shutdown(["ID1"]),(':-oncancel') = Maplets:-Elements:-
Shutdown())); DEAstep := 0; n := Dimension(f); m := ColumnDimension(g); p := 
Dimension(h); nrem := n; if not m = p then printf(
"dynext: This procedure only suitable for square systems. \n"); ERROR(
`see last comment`); end if; locxnulln := copy(xnull); xn := Vector(n,symbol =
x); y1 := Vector(m,symbol = y); locfn := copy(f); locgn := copy(g); u := Vector
(m); unull := {}; printf("Step: %a \n",0); reldegree := traperror(Reldeg(locfn,
locgn,h,x,locxnulln,u,unull,'rdeg','rtot','Adeg')); if reldegree = lasterror 
then print(reldegree); reldegree := 'false'; end if; if reldegree = true then 
printf("The relative degree of the original system well-defined. Procedure dyne\
xt is cancelled. \n"); else for i to n do printf("Step: %a \n",i); DEAstep := i
; for j to numofNonZe do if LinearAlgebra[Rank](DeleteRow(Adeg,io[j])) = 
LinearAlgebra[Rank](Adeg) then io1 := io[j]; jo1 := jo[j]; break; end if; end 
do; bdeg := h[io1]; eval(rdeg[io1]); q := 1; for k to rdeg[io1] do bdeg := 
Ldiff(locfn,bdeg,x,1); end do; nq := n+1; mq := m-1; tmpResult := false; for jj
to 10 do break if Compare("t",convert(tmpResult,string)); K := parse(Maplets:-
Display(maplet1)[]); if Compare("d",convert(K,string)) then K1 := -bdeg; else 
K1 := K; end if; S := parse(Maplets:-Display(maplet2)[]); if Compare("d",
convert(S,string)) then S1 := 1; else S1 := S; end if; xhelp := Vector(1); xn 
:= Vector(nq,symbol = x); xhelp := x[n+1]; fhelp := locfn+(`*`~(locgn[1 .. (
NULL),jo1],K1)) /~ Adeg[io1,jo1]; fnew := fhelp+((`*`~(locgn[1 .. (NULL),jo1],
S1)) *~ xhelp) /~ Adeg[io1,jo1]; for kk from n+1 to nq do printf(
"Adding state element x[%a] \n",kk); end do; fnew := Concatenate(1,fnew,Vector(
q)); g11 := IdentityMatrix(q); g12 := Matrix(n,q); gnew := Matrix(n+q,m); gnew[
1 .. (NULL),jo1] := <g12, g11>; g22 := Matrix(q,m); g21 := Matrix(n,m); for j 
to m do if j <> jo1 then g21[1 .. (NULL),j] := convert(convert(locgn[1 .. (NULL
),j],Matrix)-((convert(locgn[1 .. (NULL),jo1],Matrix) *~ Adeg[io1,j]) /~ Adeg[
io1,jo1]),Vector); end if; end do; g2 := <g21, g22>; for j to m do if j <> jo1
then gnew[1 .. (NULL),j] := g2[1 .. (NULL),j]; end if; end do; for j from n+1 
to nq do xnhelp := {x[j] = x[j]}; locxnulln := locxnulln union xnhelp; end do;
tmpU := K1+S1*xhelp; tmpU := tmpU/Adeg[io1,jo1]; for j to m do if j <> jo1 then
tmpU := tmpU-Adeg[io1,j]*u[j]; end if; end do; IntMaplet2 := Maplet([["Consider\
ing the following result, do you want to continue with chosen k(x) & s(x)"], [
"Iteration No (itr):", MathMLViewer(('value') = MathML[Export](i),('height') =
20,('width') = 50,'fontsize = 11')], ["V^(itr-1):", MathMLViewer(('value') = 
MathML[Export](simplify(tmpU)),('height') = 120,('width') = 400,'fontsize = 10'
)], ["f^(itr)=", MathMLViewer(('value') = MathML[Export](map(simplify,fnew)),('
height') = 450,('width') = 800,'fontsize = 10')], ["g^(itr)=", MathMLViewer(('
value') = MathML[Export](map(simplify,gnew)),('height') = 200,('width') = 400,'
fontsize = 10')], ["Continue with the selected functions yes/no?: ", TextField[
'TF1']()], Button("OK",Shutdown(['TF1']))]); result := Maplets[Display](
IntMaplet2); tmpResult := Compare("y",result[1]); end do; locfn := map(simplify
,convert(fnew,Vector)); locgn := map(simplify,gnew); xnew := Vector(nq,symbol =
x); reldegree := traperror(Reldeg(locfn,locgn,h,x,locxnulln,u,unull,'rdeg','
rtot','Adeg')); n := nq; if reldegree = lasterror then print(reldegree); 
reldegree := 'false'; end if; if reldegree = true then printf(
"The (vector) relative degree in step %a is well-defined. \n",i); PrintModel(x,
u,locfn,locgn); break; elif i = n then printf(); printf("   \n"); end if; end 
do; end if; if 5 < _npassed then fn := copy(locfn); end if; if 6 < _npassed 
then gn := copy(locgn); end if; if 7 < _npassed then hn := copy(h); end if; if
9 < _npassed then xnulln := copy(locxnulln); end if; if 10 < _npassed then 
Adegn := copy(Adeg); end if; if 11 < _npassed then rdegn := copy(rdeg); end if;
end proc;
