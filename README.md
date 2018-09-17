# Metric-graphs

This README is written for humans. Specifically it is written for humans who aren't great at Python/Sage. The code is very user friendly.

graph_completer.sagews is Sage code that goes along with the following mathematical paper:
"[Completing graphs to metric spaces](https://arxiv.org/abs/1706.00295)" by the 2016 Prague Ramsey Doccourse group.

You can compile it online, for example, at the free website:

    https://cocalc.com/

# Using the completer

The major function is:

    complete_partial_graph(graph, parameters)

'graph' is any graph without loops or multiedges with integer labels between 1 and delta.

'parameters' is a list of the form: 

    [delta, C, Cmax, K1, K2, ... ]

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
    [line 878] sample_graph_number = m  

Then uncomment (i.e. remove the # from) line 891.

# Show me more steps!

By default, the completer will not show steps, but you can add optional arguments to show more:

    complete_partial_graph(graph, parameters, display_all_steps, display_recursive_steps)

'display_all_steps' is True or False. It prints all of the steps in the algorithm or not. It will not show you the steps where the algorithm calls itself to complete connected components.

'display_recursive_steps' is True or False. It prints the recursive steps that 'display_all_steps' passes over. Now you see everything!

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

Either input your own parameters, or choose an with 0 < n < 19 for pre-loaded parameters. (I suggest trying n = 15.)

8.2 The choice of graph

Either input your own graph, or choose an with 0 < m < 15 for pre-loaded parameters. (I suggest trying m = 1.)

8.3 Choose whether to display all steps

The defaults are:

    display_all_steps = True
    display_recursive_steps = True

Set these to be False if you want less stuff displayed.

8.4 The function call

Uncomment this line (i.e. remove the \# ) and you're good to go!
