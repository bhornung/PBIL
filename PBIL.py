"""
Module containing a rudimentary implementation of the Population Based Incremental Learning algorithm.
(For details see: https://www.ri.cmu.edu/pub_files/pub1/baluja_shumeet_1994_2/baluja_shumeet_1994_2.pdf)
"""

import numpy as np

def tournament_selection(fitness_list, tournament_size, ntournaments):
    """
    Performs a tournament selection.
    Parameters:
        fitness_list (np.ndarray 1D) : the list of fitnesses of the individuals. It is **not** sorted!
        nchoose (int) : the tournament size. nchoose individuals will be compared at each round.
        ntournaments (int) : the number of tournaments.
    Returns:
        idcs_selected (np.ndarray if shape (nchoose)) : the indices of the selected individuals
    """
    
    idcs_selected = np.zeros(ntournaments, dtype = np.int) 
    idcs = np.arange(fitness_list.shape[0])
    
    for i in range(ntournaments):
        idcs_compare = np.random.choice(idcs, size = tournament_size, replace = False)        
        idcs_selected[i] = idcs_compare[np.argmax(fitness_list[idcs_compare])]
        
    return idcs_selected
	
	
def elitist_selection(fitness_list, n_selected):
	"""
	Chooses the fittest genes.
	Parameters:
		fitness_list (np.ndarray) : fitness list of the population
		n_selected (int) : number of individuals to be selected
	Returns:
		idcs_selected (np.ndarray of int) : indices of the fittest individuals
	"""
	return np.argsort(fitness_list)[-n_selected:]
	

class PBIL(object):
    """
    Population based incremental learning (PBIL) class.
	
	Public attributes:
		best_fitness (np.float) : best  fitness encountered during optimisation
		
		best_individual (np.ndarray) : the fittest individual encountered during optimisation
		
		ngenes (int) : number of genes
		
		optim_kwargs ({:}) : generic optimisation control parameters:
			maxiter (int ): maximum number of optimisation cycles. Default: 200.
			target_value (float) :  the expected maximum of the target function. Default: 1.0.
			learning_rate (float) : the contribution weight of the new distribution. Default:  0.02
			mutation_rate (float) : the probability of a gene being mutated. Default: 0.02
			mutation_shit (float) : the weight of the mutation
			
		popsize (int) : number of individuals
			
	Methods:
		fit() : perform at maximum maxiter optimisation cycles
    """
    
    OPTIM_KWARGS = {'maxiter' : 200, 
                    'target_value' : 0.0,
                    'learning_rate' : 0.02,
                    'mutation_rate' : 0.02,
                    'mutation_shift' : 0.02}
    
    @property
    def best_fitness(self):
        return self._best_fitness
    
    @property
    def best_individual(self):
        return self._best_individual
    
    @property
    def ngenes(self):
        return ngenes
    
    @property
    def optim_kwargs(self):
        return self._optim_kwargs
    
    @property
    def popsize(self):
        return self.popsize
    
    def __init__(self, popsize, ngenes, 
                 func, 
                 func_args = [], 
                 func_kwargs = {},
                 optim_kwargs = {}, 
                 selection_func = lambda x: x,
                 selection_args = [],
                 selection_kwargs = {}
                ):
        """
        Initialises an instance of the PBIL optimiser.
		Parameters:
			popsize (int) : number of individuals
			ngenes (int) : number of genes in an individual
			func (callable) : loss function. This function will be maximised in accordance with the GA convention.
			func_args ([]) : optional arguments of the loss function
			func_kwargs ({:}) : optional keyword arguments of the loss function
			optim_kwargs ({:}) optional keyword arguments controlling the optimisation. See class docstring for details.
			selection_func (callable) : function to select parents from the population. Default: random selection.
			selection_args ([]) : optional arguments of the selection function
			selection_kwargs ({:}) : optional keyword arguments of the selection function
        """
        
        # generic control parameters of the algorithm
        self._popsize = popsize
        self._ngenes = ngenes
        
        self._optim_kwargs = PBIL.OPTIM_KWARGS
        
        for k, v in optim_kwargs.keys():
            if k not in self._optim_kwargs.keys():
                raise KeyError("Unexpected keyword for optim_kwargs: {0}".format(k))
            else:
                self._optim_kwargs = v
        
        # loss function
        self._func = func
        self._func_args = func_args
        self._func_kwargs = func_kwargs
        
        # selection algorithm
        self._selection_func = selection_func
        self._selection_args = selection_args
        self._selection_kwargs = selection_kwargs
                   
        # optimisation outputs
        self._best_fitness = -100000000.0
        self._best_individual = None
             
        # initialise population
        self.initialise_population()
        self._distribution = self.calculate_distribution(self._population)
        
    def calculate_distribution(self, population):
        """
        Calculates the distribution of bits for each gene.
        """
        distribution = np.mean(population, axis = 0)
        
        return distribution
    
    def calculate_fitness(self):
        """
        i) Calculates fitness across the entire population.
        ii) finds best gene
        iii) stores best gene
        iv) stores best fitness
        """  
        self._fitness_list = np.apply_along_axis(self._func, 1, self._population,
                                                 *self._func_args, **self._func_kwargs)
    
        i_best = np.argmax(self._fitness_list)
        best_fitness_temp = self._fitness_list[i_best]
        
        if best_fitness_temp > self._best_fitness:
            self._best_fitness = self._fitness_list[i_best]
            self._best_individual = self._population[i_best]

    def generate_new_population(self):
        """
        Generates a new population using a binomial distribution.  
        """
        
        for idx in range(self._popsize):
            self._population[idx] = np.random.binomial(1, self._distribution)
    
    def initialise_population(self):
        """
        Initialises a new population of
        self.popsize individuals of
        self.ngenes length 
        """
        
        self._population = np.random.randint(0, high = 2, size = (self._popsize, self._ngenes))    

    def mutate_distribution(self):
        """
        Alters the population's underlying binary distribution.
        """
        
        mutation_rate = self._optim_kwargs['mutation_rate']
        mutation_shift = self._optim_kwargs['mutation_shift']
        
        # select mutation sites
        mask_mutate = np.random.choice([False, True], size = (self._distribution.shape), 
                                       p = [1.0 - mutation_rate, mutation_rate])
        
        # zero genes selected
        if not np.any(mask_mutate):
            return
        
        idcs_mutate = np.arange(mask_mutate.size)[mask_mutate]
        
        # calculate mutation
        mutation = np.random.randint(0, high = 2, size = idcs_mutate.size)
        
        self._distribution[idcs_mutate] = self._distribution[idcs_mutate] * (1.0 - mutation_shift) \
                                        + mutation * mutation_shift
        
    def refresh_distribution(self, new_distribution):
        """
        Updates the gene distribution. The new distribution is the linear combination
        of the current one and thos4e derived from the population of selected individuals.
        Updates:
            self._distribution
        """
        learning_rate = self._optim_kwargs['learning_rate']
        
        self._distribution = (1.0 - learning_rate) * self._distribution + learning_rate * new_distribution
        
    def select_parents(self):
        """
        Selects the parents from the population using
        the user defined selection algorithm.
        Returns:
            parents (np.ndarray) : a selection of individuals
        """
        
        idcs_selected = self._selection_func(self._fitness_list, 
                                             *self._selection_args,
                                             **self._selection_kwargs)
        
        parents = self._population[idcs_selected]
        
        return parents

    def fit(self):
        """
        Attempts to maximise the target function.
        """
        
        i_iter = 0
        
        while (i_iter < self._optim_kwargs['maxiter']):
                
                i_iter += 1
                
                # evaluate population
                self.calculate_fitness()
                
                print("Iteration: {0}. Fitness of best gene: {1}".format(i_iter, self._best_fitness))
                
                # shortcut on convergence
                if self._best_fitness == self._optim_kwargs['target_value']:
                    return
                
                # select best performing genes
                parents = self.select_parents()
                
                # the distribution of the genes above
                new_distribution = self.calculate_distribution(parents)
                
                # bias old distribution
                self.refresh_distribution(new_distribution)
                
                # mutate distribution
                self.mutate_distribution() 
                
                self.generate_new_population()