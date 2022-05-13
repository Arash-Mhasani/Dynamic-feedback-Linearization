

# Loading the system dynamics, and the required packages

restart:
with(LinearAlgebra):
with(ArrayTools):
with(Statistics):
with(StringTools):
with(Maplets[Elements]):
with(Matlab):
with(CodeGeneration):
interface(rtablesize=20):


# Loading the required subfunctions  

read "RelDeg.mpl":
read "dynExt.mpl":
read "printModel.mpl":
read "lDiff.mpl":
read "Models.mpl":




# Applying the dynamic extension algorithm 
# 

DynExt(f, g, h, x, X0, 'fn', 'gn', 'hn', 'xn', 'X0', 'Adegn', 'rdegn'):

NULL;




