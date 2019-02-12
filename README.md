# PBIL

Population Based Incremental Learning (PBIL) optimiser. For an enjoyable digression of PBIL pease be referred to [this article](For details see: https://www.ri.cmu.edu/pub_files/pub1/baluja_shumeet_1994_2/baluja_shumeet_1994_2.pdf). It is designed to find the optimal binary encoded solution vector.

## Class description

```python
    PBIL(popsize,
         ngenes, 
         func, 
         func_args = [], 
         func_kwargs = {},
         optim_kwargs = {}, 
         selection_func = lambda x: x,
         selection_args = [],
         selection_kwargs = {})
```

Initialises an instance of the PBIL optimiser.

### Parameters:
    
**popsize (int)** : number of individuals
      
**ngenes (int)** : number of genes in an individual
      
**func (callable)** : loss function. This function will be maximised in accordance with the GA convention.
      
**func_args ([])** : optional positional arguments of the loss function
      
**func_kwargs ({:})** : optional keyword arguments of the loss function
      
**optim_kwargs ({:})** optional keyword arguments controlling the optimisation.
      
**selection_func (callable)** : function to select parents from the population. Default: random selection.
      
**selection_args ([])** : optional arguments of the selection function
      
**selection_kwargs ({:})** : optional keyword arguments of the selection function
      
#### `optim_kwargs` 

Generic optimisation control parameters:

**maxiter (int )** : maximum number of optimisation cycles. Default: 200.
      
**target_value (float)** :  the expected maximum of the target function. Default: 1.0.
      
**learning_rate (float)** : the contribution weight of the new distribution. Default:  0.02
      
**mutation_rate (float)** : the probability of a gene being mutated. Default: 0.02
      
**mutation_shift (float)** : the weight of the mutation. Default: 0.02.
      
 ### Attributes
 
  These attributes are accessible as read only fields.
 
  **best_fitness (np.float)** : best  fitness encountered during optimisation
    
  **best_individual (np.ndarray)** : the fittest individual encountered during optimisation
    
  **ngenes (int)** : number of genes
    
  **popsize (int)** : population size
    
  **optim_kwargs** : generic optimisation control parameters
    
### Methods

#### fit()

Optimises the target function. The optimal solution vector can be retrieved by the `best_individual` attribute. 
