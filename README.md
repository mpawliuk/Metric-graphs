# Metric-graphs

graph_completer.sagews is sage code that goes along with the following mathematical paper:
"Metric graphs and their Ramsey completions" (or something like that) by the 2016 Prague Ramsey Doccourse group.

# Using the completer

The major function is:

    complete_partial_graph(graph, parameters, display_all_steps)

'graph' is any graph without loops or multiedges with integer labels between 1 and delta.

'parameters' is a list of the form: 

    [delta, C, Cmax, K1, K2, ... ]

'display_all_steps' is True or False and prints all of the steps in the algorithm or not.

# Preloaded graphs and parameters

There are 18 sets of parameters preloaded in the code (lines 738-763). To use them call:

    parameters_sample[n]

where n is between 1 and 18.

There are 14 preloaded classes of graphs in the code (lines 765-764). To use them call:

    graph_sample(n,parameters)

where n is between 1 and 14, and 'parameters' is your favourite list of parameters. Some of these graphs depend on the parameters. For example, the graph n=11 is an edge with label delta and a single point.

The preloaded parameters and graphs contain all the ones mentioned in the body of the paper.

 -----
To make you life even easier you can just choose numbers 0 < n < 15 and 0<m<19 and adjust the following lines of the code:

    [line 870] parameters_start = parameters_sample[n]
    [line 876] sample_graph_number = m  

Then uncomment (i.e. remove the # from) line 906.

# What is the rest of the code?

The code is divided into sections based on topic.

## Computer Stuff:

### 1. All about parameters

  Generate acceptable parameters, check that parameters are acceptable, 
  get info about parameters, create magic numbers
  
### 2. All about triangles

  Create all triangles, create forbidden triangles, check if graph contains a forbidden triangle
  
### 3. All about forks

  Generate forks, sort them by time, find forks in a graph
  
### 4. Antipode and bipartite stuff

  Find antipodes, make an antipodal parition, make a graph antipodally symmetric,
  how the final step works for bipartite, check to see if a graph is antipodal

### 5. Printing stuff

  Highlight edges when adding new ones, change the color of the nodes based on the bipartition
  
### 6. The major function

  The completion algorithm and printing code.
  
### 7. Testing functions
 
  Three functions designed to test many graphs and parameters at once.
 
## Human Stuff:

8.0 The list of sample graphs and parameters

8.1 The choice of parameters

8.2 The choice of graph

8.3 Choose whether to display all steps

8.4 The function call
