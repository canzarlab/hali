# Hali #
### Comparison of (rooted) phylogenetic trees and directed acyclic graphs based on arboreal matchings. ###

## Installing Hali ##

All dependencies are bundled with Hali. In particular, Hali implements its own non-linear solver and does 
not rely on any external (I)LP solver. To ```convert``` different input formats accepted by Hali you need 
to additionaly install Boost graph library. To build Hali, simply run

```
make
```

## Running Hali ##

### Usage ###

Hali can be run in three different modes. On *phylogenetic trees* in newick format:

```hali <filename.newick> <filename.newick> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>```

On *gene ontologies*, given in a format described below: 

```hali <yeastnet> <mapping> <go> <align> <constraints> <weightfunc> <k> <solver>```

On *tumor progression trees* produced by our tool convert from DOT trees. 

```hali <tree> <map> <tree> <map> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>```

### Arguments ###

Hali implements various strategies to find a (near) optimal solution. A simple greedy strategy, an exact 
branch and bound scheme

`<filename.newick>`
: input tree in newick file format

`<tree>`
  : input tree file in the following format (use convert.cpp to convert from .dot)      
  [child node] [parent node] default [newline] ...   
  

`<map>`
  : format  
  [label] [node] [newline] ...

`<yeastnet>`
  : input dag in the following format:  
  [child node] [parent node] default [newline] ...

`<mapping>`
  : input gene annotation in the following format:  
  [gene] [node] [newline] ...

`<go>`
  : input dag in the following format:  
  [parent node] [child node] default [newline] ...  
  [node] [gene] gene [newline] ...

`<align>`
  : output alignment file in the following format:  
  [node in first graph] [node in second graph] [fractional solution] [newline] ...

`<constraints>`
  : 0=matching only   
  1=forbid crossing edges  
  2=forbid crossing and semi-independent edges

`<weightfunc>`
  : j=jaccard  
  s=symdif

`<k>`
  : jaccard exponent (ignored in case of symdif)

`<vareps>`
  : weight function threshold (all edges with weight below vareps will be ignored)

`<coneps>`
  : **TODO**

`<solver>`
  : 0=greedy  
  1=geno fractional  
  2=branch and bound  
  3=covering-packing  
  4=greedy branch and bound  
  5=geno integral  
  6=geno fractional then integral  
