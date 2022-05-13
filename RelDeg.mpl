Reldeg := proc (f, g, h, x, xnull, u, unull, rdeg, rtot, Adeg, cond) local 
locrdeg, locrtot, lgh, lfh, setvar, i, j, k, ii, jj, kk, varseq, solseq, 
notdefseq, def, hulpvar, lghaug, condhlp, condseq, locxnull, temp, locunull, 
xwork, rddef, Adegi, value, notdef, notaff, printbd, printbd2, printbd3, 
printrdeg, reldeg, zinput, zc, n, p, Adegnull, mnull, io1, jo1; global m, io, 
jo, numofNonZe, DEAstep; if _npassed < 5 and 12 < _npassed then ERROR(
`reldeg: invalid number of arguments`); end if; n := Dimension(f); p := 
Dimension(h); m := p; locxnull := xnull; locunull := unull; xwork := locxnull;
if not p <= m then printf(
"reldeg: number of inputs must be larger than or equal to number of outputs. \n\
"); printf("   \n"); end if; locrdeg := Vector(p); lfh := Vector(p); lgh := 
Matrix(1 .. p,1 .. m); varseq := seq(x[i],i = 1 .. n); hulpvar := n; locrtot :=
0; rddef := 'true'; notaff := NULL; notdef := NULL; for i to p do notdefseq :=
NULL; def := 'true'; condseq[i] := NULL; setvar := 0; value := n; for k to 
value while setvar = 0 do locrdeg[i] := k; if k = 1 then lfh[i] := h[i]; else 
lfh[i] := Ldiff(f,lfh[i],x,1); end if; for j to m do lgh[i,j] := simplify(
normal(Ldiff(g[1 .. (NULL),j],lfh[i],x,1))); end do; for j to m do if lgh[i,j]
<> 0 then setvar := 1; traperror(simplify(normal(eval(subs(xwork,lgh[i,j])))));
if % = lasterror then printf("reldeg: WARNING current working point can not be \
substituted in the (internal) computations of the relative degree but the compu\
tations are proceeded until a fatal ERROR will occur \n "); elif % = 0 then def
:= 'false'; notdefseq := notdefseq, j; else def := 'true'; break; end if; end 
if; end do; if not def then notdef := notdef, i; notdefseq := [notdefseq]; for
kk to nops(notdefseq) do solseq := lasterror; if solseq <> lasterror then for 
jj to nops(solseq) do condhlp := NULL; for kk to nops(op(jj,solseq)) do if not
op(1,op(kk,op(jj,solseq))) = op(2,op(kk,op(jj,solseq))) then condhlp := condhlp
, op(kk,op(jj,solseq)); end if; end do; if condhlp = NULL then condhlp := all;
end if; condseq[i] := condseq[i], condhlp; end do; end if; end do; condseq[i] 
:= {condseq[i]}; end if; if k = value and setvar = 0 then rddef := 'false'; 
notaff := notaff, i; end if; end do; locrtot := locrtot+locrdeg[i]; end do; 
printrdeg := NULL; reldeg := NULL; if 8 < _npassed then rtot := locrtot; end if
; if 9 < _npassed then rdeg := locrdeg; Adeg := lgh; printf(
"The decoupling matrix is: \n"); print(A^cat(`<`,DEAstep,`>`) = eval(Adeg)); 
printf("      \n"); printf("      \n"); printf("      \n"); end if; if 0 < nops
({notaff}) then rddef := 'false'; printf("reldeg: The relative degree cannot be\
 defined for this system, the output(s) below are not at all affected by the in\
put but are only depending on the initial state. \n "); printf(
"reldeg: The rank of the matrix Adeg is less than the number of output(s). \n")
; print(not_affected_outputs = eval(notaff)); printf("      \n"); printf(
"      \n"); printf("      \n"); for i to p do if has({notaff},i) then rdeg[i]
:= nd; end if; end do; rtot := nd; printrdeg := not_defined; reldeg := 
not_defined; end if; if 0 < nops({notdef}) then rddef := 'false'; Adegnull := 
map(simplify,map(normal,subs(xwork,op(Adeg)))); numofNonZe := 0; io := Vector(m
); jo := Vector(m); for j to m do for i to m do if Adegnull[i,j] <> 0 then 
numofNonZe := numofNonZe+1; io[numofNonZe] := i; jo[numofNonZe] := j; end if; 
end do; end do; for j to numofNonZe do if LinearAlgebra[Rank](DeleteRow(Adeg,io
[j])) = LinearAlgebra[Rank](Adeg) then io1 := io[j]; jo1 := jo[j]; break; end 
if; end do; printf("The evaluation of the decoupling matrix at x0 is: \n"); 
print(nprintf(`#msup(mi("A"),mn("<%d>"))`,DEAstep)(nprintf(
`#msup(mi("x0"),mn("<%d>"))`,DEAstep)) = Adegnull); mnull := LinearAlgebra:-
Rank(Adegnull); printf("The rank of the above matrix at x0 is less than the num\
ber of outputs, %a. Therefore, the relaive degree is not defined at x0. \n",m);
printf("We select indices: \n"); print((i__0, j__0) = (io1, jo1)); printf(
"      \n"); printf("      \n"); printf("      \n"); printrdeg := not_defined;
reldeg := not_defined; end if; if Determinant(map(simplify,map(normal,subs(
xnull,Adeg)))) = 0 and nops({notdef}) = 0 and nops({notaff}) = 0 then rddef :=
'false'; printf("reldeg: Although each output has a well defined relative degre\
e the  of the matrix Adeg is less than the number of output(s) due to the fact \
that the following input(s) appear always later than other(s): \n "); zinput :=
NULL; for i to m do zc[i] := NULL; for j to p do zc[i] := zc[i], lgh[j,i]; end
do; if {zc[i]} = {0} then zinput := zinput, i; end if; end do; print(
late_inputs = zinput); printf(`     `); printrdeg := not_defined; reldeg := 
not_defined; Adegnull := map(simplify,map(normal,subs(xwork,op(Adeg)))); 
numofNonZe := 0; io := Vector(m); jo := Vector(m); for j to m do for i to m do
if Adegnull[i,j] <> 0 then numofNonZe := numofNonZe+1; io[numofNonZe] := i; jo[
numofNonZe] := j; end if; end do; end do; elif lgh = n and locrtot = n and 
reldeg = NULL then printf("The total relative degree at the working point is eq\
ual to the number of states. \n"); end if; printf("r-indices are given by: ");
print(r^cat(`<`,DEAstep,`>`) = eval(locrdeg^%T)); if nops({notdef}) <> 0 or 
nops({notaff}) <> 0 then printf("      \n"); printf("      \n"); printf(
"      \n"); end if; if Determinant(map(simplify,map(normal,subs(xnull,lgh))))
= 0 and nops({notdef}) = 0 and nops({notaff}) = 0 then printf("The number(s) be\
tween [] give the vector and total relative degree of each seperate output, not\
_defined means that this relative degree is not properly defined. \n"); printf(
"      \n"); printf("      \n"); printf("      \n"); end if; rddef; end proc;
