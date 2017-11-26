# Hali #
### Comparison of (rooted) phylogenetic trees and directed acyclic graphs based on arboreal matchings. ###

## Installing Hali ##

All dependencies are bundled with Hali. In particular, Hali implements its own non-linear solver and does 
not rely on any external (I)LP solver. To ```convert``` different input formats accepted by Hali you need 
to install the Boost graph library. To build Hali, simply run

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

Hali implements various strategies to find a (near) optimal solution. Its non-linear solver is based on an 
augmented Lagrangian approach. Hali can provide an optimal fractional solution or an optimal integral solution
obtained through branch and bound. It can also find a (suboptimal) integral solution based on a greedy 
strategy or enforce integrality by non-linear constraints. 

`<solver>`
  : 0=greedy  
  1=fractional  
  2=branch and bound  
  3=covering-packing  
  4=greedy branch and bound  
  5=non-linear integral  
  6=warm-start integral from fractional  

Options 3,4 and 6 combine different solver strategies that require problem specific tuning.

#### Input/Output formats        
     
`<filename.newick>`
  : input tree in newick file format

`<tree>`
  : input tree file in the following format (use convert.cpp to convert from .dot):        
  [child node] [parent node] default [newline] ...   
  
`<map>`
  : format [label] [node] [newline] ...  
  A label can encode, for example, the copy number of genes probed in FISH data from single cells.

`<yeastnet>`
  : input DAG in the following format:  
  [child node] [parent node] default [newline] ...

`<mapping>`
  : input gene assignment in the following format:  
  [gene] [node] [newline] ...

`<go>`
  : input DAG (GO) in the following format:  
  [parent node] [child node] default [newline] ...  
  [node] [gene] gene [newline] ...

`<align>`
  : output the final alignment in file <align> in the following format:  
  [node in first graph] [node in second graph] [fractional solution] [newline] ...
  
#### Model adjustment      
        
`<constraints>`
  : 0=unconstrained matching   
  1=forbid crossing edges  
  2=forbid crossing and semi-independent edges

`<weightfunc>`
  : j=jaccard  
  s=symmetric difference

`<k>`
  : order of jaccard weight (ignored in case of symmetric difference )

#### Miscellaneous   
   
`<vareps>`
  : edge weight threshold (all edges with weight below vareps will be ignored) - default 0

`<coneps>`
  : **TODO**
