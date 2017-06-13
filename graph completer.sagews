###############################################################################
### Stuff for computers  ######################################################
###############################################################################

#########################
### 1. For parameters ###
#########################

def get_C(parameters):
    """Get C0 and C1.

    C0 is even, C1 is odd.
    """
    C, Cmax = parameters[1:3]
    C_list = [0,0]

    if C % 2 == 0:
        C_list[0] = C
        C_list[1] = Cmax
    else:
        C_list[0] = Cmax
        C_list[1] = C
    return C_list

def is_acceptable(parameters):
    """ Check if parameters are acceptable.

    See Definition 3.1.
    """
    delta, C, Cmax, K1, K2 = parameters[0:5]
    C0 = get_C(parameters)[0]
    C1 = get_C(parameters)[1]
    if delta < 2:
        return False
    if (K1 > 3*delta):
        if (K2 > 0) or (C1 != 2*delta + 1):
            return False
        else:
            return True
    if not (K1 in range(1,2*delta + 1)):
        return False
    if not (K2 in range(1,2*delta + 1)):
        return False
    if not (C0 in range(2*delta + 1, 3*delta+3)):
        return False
    if not (C1 in range(2*delta + 1, 3*delta+3)):
        return False
    return True

def what_case(parameters):
    """Return a dictionary that captures all info about the case.

    --"case" in [1,2a,2b,3]
    --"antipodal" in [True, False]
    --"bipartite" in [True, False]
    --"even" in [True, False]
    --"odd" in [True, False]
    --"subcase" in ["3 constrained", "Bipartite 3 contrained", 1,2,3,4]
    ----1. antipodal classes in Case (IIA) with even d (Corollary 7.5),
    ----2. antipodal bipartite classes in Case (I) with odd d (Corollary 7.7),
    ----3. antipodal classes in Case (IIA) with odd d (Theorem 7.4), and
    ----4. antipodal bipartite classes in Case (I) with even d (Theorem 7.6).
    """
    delta, C, Cmax, K1, K2 = parameters[0:5]
    delta_parity = delta % 2

    if not is_acceptable(parameters):
        print "These parameters are not acceptable.", parameters
        return None

    # Case 1 - Bipartite
      # K1 > 3*delta is used instead of infinity.
    if (K1 > 3*delta):
        if (Cmax > 2*delta + 3):
            return {'case': '1', 'antipodal':False, 'bipartite':True,
                    'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                    'subcase': 'Bipartite 3 constrained'}
        elif (Cmax == 2*delta + 2):
            if (delta_parity == 1):
                return {'case': '1', 'antipodal':True, 'bipartite':True,
                        'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                        'subcase': 2}
            else:
                return {'case': '1', 'antipodal':True, 'bipartite':True,
                        'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                        'subcase': 4}

    # Case 2
    if (C == 2*K1 + 2*K2 + 1) and (K1 + K2 >= delta) and (K1 + 2*K2 <= 2*delta - 1):
        # Case 2a
        if Cmax == C + 1:
            # Case 2a Antipodal
            if (C == 2*delta + 1):
                if (delta_parity == 0):
                    return {'case': '2a', 'antipodal':True, 'bipartite':False,
                            'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                            'subcase': 1}
                else:
                    return {'case': '2a', 'antipodal':True, 'bipartite':False,
                            'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                            'subcase': 3}
            else:
                return {'case': '2a', 'antipodal':False, 'bipartite':False,
                        'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                        'subcase': '3 constrained'}
        # Case 2b
        elif (K1 == K2) and (3*K2 == 2*delta - 1):
            return {'case': '2b', 'antipodal':False, 'bipartite':False,
                    'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                    'subcase': '3 constrained'}

    # Case 3
    if (C > 2*delta + K1) and (K1 + 2*K2 >= 2*delta - 1) and (3*K2 >= 2*delta):
        if (K1 + 2*K2 == 2*delta - 1) and (C < 2*delta + K1 + 2):
            return "Case 3, Special Case 1 Violated"
        if (Cmax > C + 1) and (C < 2*delta + K2):
            return "Case 3, Special Case 2 Violated"
        else:
            return {'case': '3', 'antipodal':False, 'bipartite':False,
                    'even': (delta_parity == 0), 'odd':(delta_parity == 1),
                    'subcase': '3 constrained'}

    print "These parameters are not admissible."
    return None

def min_magic_number(parameters):
    """Return the minimal magic number if it exists."""
    delta, C, Cmax, K1, K2 = parameters[0:5]
    delta_parity = parameters[0]%2
    case_info = what_case(parameters)

    # x // 2 takes floor of divsion.
    # -(-x // 2) takes ceiling of division.
    if case_info['antipodal']:
        if delta_parity == 0:
            M1 = M2 = delta // 2
        elif delta_parity == 1:
            M1 = M2 = (delta - 1) // 2
    elif case_info['bipartite']:
        M1 = delta // 2
        M2 = (Cmax - delta - 1) // 2
    elif case_info['subcase'] == "3 constrained":
        M1 = max(K1, -(-delta // 2))
        M2 = min(K2, (C - delta - 1) // 2)
    else:
        # This means that something has gone wrong.
        M1, M2 = 1, 0

    assert M1 <= M2
    return M1

##############################
### 2. All about triangles ###
##############################

def make_triangles(delta):
    """Return all possible triangles (good and bad).

    Output format is [length 1, length 2, length 3, "Good/bad"]
    Written in nondecreasing order.
    """
    triangles = []
    for i in range(1,delta+1):
        for j in range(i,delta+1):
            for k in range(j,delta+1):
                triangles.append([i,j,k,""])
    return triangles

def make_forbidden_triangles(parameters):
    """Return all forbidden triangles."""
    delta, C, Cmax, K1, K2 = parameters[0:5]
    all_triangles = [x for x in make_triangles(delta)]

    C0 = get_C(parameters)[0]
    C1 = get_C(parameters)[1]

    forbidden_triangles = [] # This is where we store the bad triangles.

    for triangle in all_triangles:
        p = triangle[0]+triangle[1]+triangle[2] # The perimeter of the triangle.
        p_parity = p % 2 # The parity of the perimeter.
        m = triangle[0] # This is the minimum side length, it's already sorted!

        # Check the triangle inequality.
        # (only need to check that the biggest side can't be reached)
        if (triangle[2] > triangle[1]+triangle[0]): #Triangle inequality
            triangle[3] = "Metric"
            forbidden_triangles.append(triangle)
        elif (p_parity == 0) and not (p < C0): # Even condition
            triangle[3] = "Even"
            forbidden_triangles.append(triangle)
        elif (p_parity == 1):
            if not (p < C1):
                triangle[3] = "Odd, C1"
                forbidden_triangles.append(triangle)
            elif not (2*K1 + 1 <= p):
                triangle[3] = "Odd, K1"
                forbidden_triangles.append(triangle)
            elif not (p < 2*K2 + 2*m):
                triangle[3] = "Odd, K2"
                forbidden_triangles.append(triangle)
    return forbidden_triangles

def contains_forbidden_triangles(graph, parameters):
    """Check if the graph has a forbidden triangle using brute force."""
    delta, C, Cmax, K1, K2 = parameters[0:5]
    forb = [x[:3] for x in make_forbidden_triangles(parameters)]
    vertices = range(len(graph.vertices()))
    for i in vertices:
        for j in vertices[i+1:]:
            for k in vertices[j+1:]:
                vi,vj,vk = graph.vertices()[i], graph.vertices()[j], graph.vertices()[k]
                lab1 = int(graph.edge_label(vi,vj))
                lab2 = int(graph.edge_label(vj,vk))
                lab3 = int(graph.edge_label(vi,vk))
                if sorted([lab1, lab2, lab3]) in forb:
                    return ((i,j,k), sorted([lab1, lab2, lab3]))
    return False

def is_complete(graph):
    """Check if a labelled grpah is complete."""
    return not graph.complement().edges(labels=False)

####################
### 3. For forks ###
####################

def generate_forks(parameters):
    """Return a list with elements of the form [(i,j),a,b,c, *_ ].

    a,b,c,*_ are the ways to complete the fork i,j.
    """
    forb = [x[:3] for x in make_forbidden_triangles(parameters)]
    delta = parameters[0]
    forks = []

    for i in range(1,delta+1):
        for j in range(i,delta+1):
            acceptable=[(i,j)]
            for k in range(1,delta+1):
                if sorted([i,j,k]) not in forb:
                    acceptable.append(k)
            if len(acceptable)>1:
                forks.append(acceptable)
    return forks

def sort_forks(parameters):
    """Return a list of the form [F(0),F(1), ..., F(delta)].

    -- F(i) will be the list of things to be completed by distance i.
    See paragraph after Definition 4.2
    """
    delta = parameters[0]
    M = min_magic_number(parameters)
    forks_by_completion = [[] for i in range(delta+1)]

    for fork in generate_forks(parameters):
        l = fork[0][0] #left fork
        r = fork[0][1] #right fork, note r>=l

        #Is the magic number in there?
        if M in fork[1:]:
            forks_by_completion[M].append(fork[0]) # Magic distance
        elif l+r < M and l+r in fork[1:]:
            forks_by_completion[l+r].append(fork[0]) # Small geodesic
        elif r-l > M and r-l in fork[1:]:
            forks_by_completion[r-l].append(fork[0]) # Big geodesic
        else:
            forks_by_completion[max(fork[1:])].append(fork[0]) # C bound, add maximal thing allowed
    return forks_by_completion

def find_forks(graph, i,j):
    """Find all forks in graph with labels i,j.

    The output elements are of the form [v,a,b],
    where v is the centre of the fork and there is no edge between a and b
    """
    forks = []
    for edge in graph.edges():
        #First find an edge with label i
        if edge[2] == str(i):
            for endpoint in edge[:2]:
                #Now look for an edge with label j coming from one of the endpoints
                for nghbr in graph.neighbors(endpoint):
                    if graph.has_edge(endpoint,nghbr, str(j)):
                        P = graph.subgraph(list(edge[:2])+[nghbr])
                        if len(P.edges()) == 2:
                            #now keep track of the missing edge
                            for vP in P.vertices():
                                if len(P.neighbors(vP)) == 2:
                                    info = [vP] + P.neighbors(vP)
                                    if info not in forks:
                                        forks.append(info)
    for edge in graph.edges():
        #Find an edge with label j
        if edge[2] == str(j):
            for endpoint in edge[:2]:
                #Now look for an edge with label j coming from one of the endpoints
                for nghbr in graph.neighbors(endpoint):
                    if graph.has_edge(endpoint,nghbr, str(i)):
                        P = graph.subgraph(list(edge[:2])+[nghbr])
                        if len(P.edges()) == 2:
                            #now keep track of the missing edge
                            for vP in P.vertices():
                                if len(P.neighbors(vP)) == 2:
                                    info = [vP] + P.neighbors(vP)
                                    if info not in forks:
                                        forks.append(info)
    return forks

def time_sort(parameters):
    """Return a list of the order in which we add edges.

    See comments after Definition 4.3 and see Figure 3.
    See also Definitions 6.2, 7.5, 7.6
    """
    delta = parameters[0]
    M = min_magic_number(parameters)
    case_info = what_case(parameters)

    temp_list = []
    for x in range(1,M):
        temp_list.append([x, 2*x - 1])
    for x in range(M+1, delta):
        temp_list.append([x, 2*(delta - x)])

    # In some cases we don't want M+1 or M-1
    if (case_info['subcase'] == 3):
        x = M + 1
        if [x, 2*(delta - x)] in temp_list:
            temp_list.remove([x, 2*(delta - x)])
    elif case_info['bipartite']:
        x = M + 1
        if [x, 2*(delta - x)] in temp_list:
            temp_list.remove([x, 2*(delta - x)])
        if case_info['antipodal'] and case_info['even']:
            x = M - 1
            if [x, 2*x - 1] in temp_list:
                temp_list.remove([x, 2*x - 1])
    return [i[0] for i in sorted(temp_list, key=lambda y:y[1])]

#######################################
### 4. Antipode and bipartite stuff ###
#######################################

def is_delta_graph(graph, parameters):
    """Check if a graph has no labels greater than delta."""
    if graph.edge_labels():
        return max([int(x) for x in graph.edge_labels()]) <= parameters[0]
    else:
        return True

def is_antipodal(graph, parameters):
    """Check whether a graph could be completed to an antipodal graph.

    Check whether any node has two different neighbours at distance delta.
    """
    delta = str(parameters[0])
    for v in graph.vertices():
        #Count the number of neighbors of distance delta
        if len([w for w in graph.neighbors(v) if graph.edge_label(v,w) == delta]) > 1:
            return False
    return True

def antipode(graph, parameters, v_input):
    """Return a vertex's antipode.

    If the vertex has multiple neighbours of max distance, it just gives one.
    """
    delta = str(parameters[0])
    for w in graph.neighbors(v_input):
        if graph.edge_label(v_input,w) == delta:
            return True, w
    return False, None

def antipodal_partition(antipodal_graph, parameters):
    """Make a choice of one vertex from each antipodal pair."""
    P = []
    for v in antipodal_graph.vertices():
        if (not (v in P)) and not (antipode(antipodal_graph, parameters, v)[1] in P):
            P.append(v)
    return P

def make_antipodally_symmetric(G_input, parameters):
    """Return a graph with two properties.

    1. Each node has an antipode.
    2. Each pair of disjoint edges of length delta are contined in a symmetric rectangle.
    See Defn 7.4.
    """
    G_output = copy(G_input)
    delta = parameters[0]
    if not is_antipodal(G_input, parameters):
        print "This graph has a vertex with two neighbours at antipodal distance. Things will not work properly."
        return Graph()
    # 1. Give each node that needs it an antipode.
    for v in G_input.vertices():
        if not antipode(G_input, parameters, v)[0]:
            #Add a vertex of distance delta
            G_output.add_edge(v, max(G_output.vertices())+1, str(delta))
    # 2. Make symmetric all antipodal quadruples. Complete all rectangles
    for edge_1 in G_output.edges():
        for edge_2 in G_output.edges():
            if (edge_1 != edge_2) and (edge_1[2] == str(delta)) and (edge_2[2] == str(delta)):
                #Make symmetric.
                for i in range(2):
                    for j in range(2):
                        v,a = edge_1[i], edge_1[1-i]
                        w,b = edge_2[j], edge_2[1-j]
                        if G_output.has_edge(v,w) and not G_output.has_edge(a,b):
                            G_output.add_edge(a,b, G_output.edge_label(v,w))
                #Complete, if a rectangle.
                for i in range(2):
                    for j in range(2):
                        v,a = edge_1[i], edge_1[1-i]
                        w,b = edge_2[j], edge_2[1-j]
                        if G_output.has_edge(v,w) and not (G_output.has_edge(v,b) or G_output.has_edge(a,w)):
                            G_output.add_edge(v, b, str(delta-int(G_output.edge_label(v,w))))
                            G_output.add_edge(a, w, str(delta-int(G_output.edge_label(v,w))))
    return G_output

def bipartite_connect(graph, parameters):
    """Connect a bipartite graph as in Lemma 6.2."""
    case_info = what_case(parameters)
    M = min_magic_number(parameters)
    components = graph.connected_components_subgraphs()

    for i, comp1 in enumerate(components):
        for comp2 in components[i+1:]:
            # Choose representatives (this is non-canonical).
            min1 = min(comp1.vertices())
            min2 = min(comp2.vertices())
            # Something special for subcase 4
            if False and (case_info['subcase'] == 4 and comp2 == components[i+1]):
                edge_parity = 0
                edge_label = str(M + edge_parity % 2)
                graph.add_edge((min1, min2, edge_label))
            else:
                for v1 in comp1:
                    for v2 in comp2:
                        edge_parity = M % 2
                        if min1 != v1:
                            edge_parity += int(comp1.edge_label(min1,v1))
                        if min2 != v2:
                            edge_parity += int(comp2.edge_label(min2,v2))
                        edge_label = str(M + edge_parity % 2)
                        graph.add_edge((v1, v2, edge_label))
    return graph

def bipartite_final_step(graph, parameters,x,y):
    """Add M or M+1 depending on if d+(x,y) = M mod 2."""
    M = min_magic_number(parameters)
    dplus = graph.shortest_path_length(x,y,by_weight=True, weight_function=lambda edge: int(edge[2]))

    # Use definition 6.2
    fix = (M+dplus)%2 #You add M or M+1
    graph.add_edge((x,y,str(M+fix)))
    return graph

#######################
### 5. For Printing ###
#######################

def vertex_colour_dict(graph, parameters):
    """Colour vertices if bipartite."""
    vertex_colour_dict = {}
    if what_case(parameters)['bipartite'] and graph.vertices():
        x0 = min(graph.vertices())
        vertex_colour_dict[rainbow(2)[0]] = [x0]
        vertex_colour_dict[rainbow(2)[1]] = []
        for v in graph.vertices():
            if v != x0:
                vertex_colour_dict[rainbow(2)[(int(graph.edge_label(x0,v))%2)]].append(v)
    return vertex_colour_dict

def display_highlighted_edges(graph, delta, edge_label):
    """Show a graph with newly added edges in black."""
    d = {str(i):rainbow(delta + 1)[i] for i in range(1, delta + 1)}
    d[edge_label] = 'black'
    graph.graphplot(edge_labels=True, color_by_label=d, layout='circular').show()
    return None

#############################
### 6. The major function ###
#############################

def complete_partial_graph(G_input, parameters, display_all_steps=False, display_recursive_steps=False):
    """Return the completed graph, subject to the parameters.

    1. Deal with Antipodal stuff.
    --1.1. Make antipodally symmetric if needed.
    --1.2. Complete one pode with reduced parameters.
    --1.3. Make antipodally symmetric if needed.
    2. Deal with bipartite stuff.
    --2.1. Complete individual components.
    --2.2. Connect components.
    3. Complete forks as perscribed.
    4. Add magic distance (+- 1) to finish.
    --4.1. Print step 4.
    """
    G_output = copy(G_input)

    delta, C, Cmax, K1, K2 = parameters[0:5]
    delta_parity = parameters[0]%2
    case_info = what_case(parameters)
    M = min_magic_number(parameters)

    if display_all_steps and not is_delta_graph(G_input, parameters):
        print "The input graph has labels larger than ", delta, ". It will not run."
        return Graph()

    ####################################
    ### 1.1. Antipode closure.       ###
    ### 1.2. Complete one pode.      ###
    ### 1.3. Antipode closure again. ###
    ####################################
    if case_info['antipodal']:
        # 1.1 Make the graph antipodally symmetric as in Defn 7.4
        G_output = make_antipodally_symmetric(G_output, parameters)

        if display_all_steps and G_output != G_input:
            print "1.1 Make graph antipodally symmetric:"
            display_highlighted_edges(G_output, delta, str(delta))
        elif display_all_steps:
            print "Graph is already antipodally symmetric."

        # 1.2. Complete one of the antipode parts as in Defn 7.5
        P = antipodal_partition(G_output, parameters)
        partition_parameters = copy(parameters)
        # Reduce delta.
        partition_parameters[0] = parameters[0] - 1
        if case_info['subcase'] in (2,4):
            # If also bipartite, reduce C to stay bipartite
            partition_parameters[1] = parameters[1] - 2
        # Add the completion of one of the parts.
        if display_all_steps:
            print "1.2 Complete the pode", P
        G_output = G_output.union(complete_partial_graph(G_output.subgraph(P), partition_parameters, \
                                                         display_all_steps, display_recursive_steps))
        if display_all_steps:
            display_highlighted_edges(G_output, delta, str(delta))

        # 1.3. Make it antipodally symmetric
        G_output = make_antipodally_symmetric(G_output, parameters)
        if display_all_steps:
            print "1.3 Make graph antipodally symmetric:"
            display_highlighted_edges(G_output, delta, str(delta))
    ##################################################
    ### 2.1. Complete each component of bipartite. ###
    ### 2.2. Bipartite connect.                    ###
    ### 2.3. Antipode closure again (if needed)    ###
    ##################################################
    if case_info['bipartite'] and (not G_output.is_connected()):
        G_output_real = copy(G_output)

        # 2.1. Complete the individual components
        if display_all_steps:
            print "2.1 Individual components will be completed."
        for comp in G_output.connected_components_subgraphs():
            G_output_real = G_output_real.union(complete_partial_graph(comp, parameters, display_recursive_steps))
        if display_all_steps:
            display_highlighted_edges(G_output_real, delta, str(delta))

        # 2.2. Now connect the parts as in Lemma 6.2
        G_output_real = bipartite_connect(G_output_real, parameters)
        if display_all_steps:
            print "2.2 Disconnected components were connected."
            display_highlighted_edges(G_output_real, delta, str(delta))

        # 2.3. Antipodally close again.
        if case_info['subcase'] == 4:
            G_output_real = make_antipodally_symmetric(G_output_real,parameters)
            if display_all_steps:
                print "2.3 Disconnected components were connected."
                display_highlighted_edges(G_output_real, delta, str(delta))
        if is_complete(G_output_real):
            return G_output_real
        else:
            # This is part of a subcase 4 run.
            # Once it is antipodally closed it will be complete.
            G_output = G_output_real

    ######################################
    ### The main part of the algorithm ###
    ######################################

    #########################
    ### 3. Complete forks ###
    #########################

    # This needs to be ordered by time
    for t in time_sort(parameters):
        added_something = False
        for fork_rule in sort_forks(parameters)[t]:
            # Things that get value t
            for fork_to_fill in find_forks(G_output,fork_rule[0],fork_rule[1]):
                # fork_to_fill looks like [v, a, b] where v is the vertex of the fork
                G_output.add_edge((fork_to_fill[1],fork_to_fill[2],str(t)))
                added_something = True
        if display_all_steps and added_something:
            print "Added edges with weight: ", t
            display_highlighted_edges(G_output, delta, str(t))
        elif display_all_steps:
            print "No edges added with weight: ", t
    #######################################
    ### 4. The final step adds M or M+1 ###
    #######################################
    for pair_of_vertices in G_output.complement().edges(labels=False):
        x = pair_of_vertices[0]
        y = pair_of_vertices[1]
        if (case_info['subcase'] in ('3 constrained', 1)):
            # For every other non-edge add the Magic distance
            G_output.add_edge((x,y,str(M)))
        else:
            # Must be bipartite case
            G_output = bipartite_final_step(G_output, parameters,x,y)
    ##################################
    #### 4.1 Print steps if needed ###
    ##################################
    if display_all_steps and not is_complete(G_output):
        if (case_info['subcase'] in ('3 constrained', 1)):
            print "Add the magic distances " +str(M)+ " to finish."
        elif (case_info['subcase'] in ('Bipartite 3 constrained', 2,3)):
            print "Add the magic distances " +str(M)+ " and " + str(M+1) + " to finish."
    ##############################
    ### Output completed graph ###
    ##############################
    return G_output

def display_graph_completion(G_input, parameters, display_all_steps, display_recursive_steps):
    """Print the run of the algorithm."""
    if display_all_steps:
        print "The parameters are:"
        print "   delta = ", parameters[0]
        print "   C = ", parameters[1]
        print "   C' = ", parameters[2]
        print "   K1 = ", parameters[3]
        print "   K2 = ", parameters[4]
        print "--------------------"
        print "The paramaters are in case: ", what_case(parameters)['case']
        print "What subcase?", what_case(parameters)['subcase']
        print "bipartite?", what_case(parameters)['bipartite']
        print "antipodal?", what_case(parameters)['antipodal']
        print "Is the diamter delta even?", what_case(parameters)['even']
        print "Is the diamter delta odd?", what_case(parameters)['odd']
        print "--------------------"
        print "Is the starting graph connected?", G_input.is_connected()
        print "--------------------"
        if what_case(parameters)['subcase'] == '3 constrained' or what_case(parameters)['subcase'] == 'Bipartite 3 constrained':
            print "Edges of a given distance will be added in order: ", time_sort(parameters)
            print "then the magic distance: ", min_magic_number(parameters), "(possibly + 1)."
            print "--------------------"

    print "Here is the starting graph:"
    colours = {str(i):rainbow(parameters[0]+1)[i] for i in range(1,parameters[0]+1)}
    G_input.graphplot(edge_labels=True, color_by_label=colours, layout='circular').show()

    G_complete = complete_partial_graph(G_input, parameters, display_all_steps, display_recursive_steps)

    if what_case(parameters)['antipodal'] and display_all_steps:
        P = antipodal_partition(G_complete, parameters)
        print "The vertices were antipodal partitioned into:", P
    if contains_forbidden_triangles(G_complete, parameters):
        print "Uh oh! The completed graph has a forbidden triangle!"
        print "The vertices: ", contains_forbidden_triangles(G_complete, parameters)[0]
        print "make the forbidden triangle: ", contains_forbidden_triangles(G_complete, parameters)[1]
    else:
        print "This graph has no forbidden triangles."

    print "--------------------"
    print "Here is the completed graph:"
    return G_complete.graphplot(edge_labels=True, color_by_label=colours, \
                                vertex_colors=vertex_colour_dict(G_complete, parameters), \
                                layout='circular').show()

############################
### 7. Testing functions ###
############################

def test_the_graphs(testing_graphs,testing_parameters):
    """Check if some subset of the sample graphs can be completed using some subset of the sample parameters.

    testing_graphs is a sublist of range(1,15).
    testing_parameters is a sublist of range(1,19).
    Especially handy for antipodal parameters.
    """
    for j in testing_parameters:
        for i in testing_graphs:
            parameters = parameters_sample[j]
            H = graph_sample(i, parameters)
            if is_delta_graph(H,parameters):
                G = complete_partial_graph(H, parameters)
                if contains_forbidden_triangles(G, parameters):
                    print i, parameters, contains_forbidden_triangles(G, parameters)
    return "All graphs have been checked."

def test_the_graphs_2(testing_parameters, extra_vertex=False):
    """Check if all paths of length 2 can be completed using some subset of the sample parameters.

    testing_parameters is a sublist of range(1,19)
    Especially handy for antipodal parameters.
    """
    H = Graph()
    H.add_edge((0,1,"1"))
    H.add_edge((1,2,"1"))
    if extra_vertex:
        H.add_vertex(3)
    for k in testing_parameters:
        parameters = parameters_sample[k]
        delta = parameters[0]
        for i in range(1,delta+1):
            for j in range(i, delta+1):
                H.set_edge_label(0,1,str(i))
                H.set_edge_label(1,2,str(j))
                if is_delta_graph(H,parameters):
                    G = complete_partial_graph(H, parameters)
                    if contains_forbidden_triangles(G, parameters):
                        print i, parameters, contains_forbidden_triangles(G, parameters)
    return "All paths have been checked."

def display_the_graphs(which_graphs, parameters):
    """Show a subset of the sample graphs.

    testing_graphs is a sublist of range(1,15).
    Note that many of the graphs depend on delta.
    """
    for j in testing_graphs:
        G = graph_sample(j,parameters)
        print j
        colours = {str(i):rainbow(parameters[0]+1)[i] for i in range(1,parameters[0]+1)}
        G.graphplot(edge_labels=True, color_by_label=colours, layout='circular').show()
    return "All graphs have been displayed."

############################################################################
### 8. Stuff for humans  ###################################################
############################################################################

parameters_sample = [
    # [Diameter, C, C', K1, K2, "Case name"]
    [5, 13, 16,   3, 3, "Case 2b"], # Smallest 2b, Example 1,2, Figure 6
    [5, 13, 16,   3, 3, "Case 2b"], # Smallest 2b, Example 1,2, Figure 6
    [8, 21, 24,   5, 5, "Case 2b"], # next smallest 2b
    [5, 13, 14,   3, 3, "Case 2a"],
    [5, 12, 13,   1, 5, "Case 3"], # No large (p>=14) cycles

    [5, 16, 17,   1, 5, "Case 3"], # Only metric forbidden
    [5, 16, 17,   5, 5, "Case 3"], # No small (p<=9) odd cycles
    [3,  7,  8,   1, 2, "Case 2a - Antipodal, non-bipartite, odd delta"], # Figure 10,
    [4,  9, 10,   1, 3, "Case 2a - Antipodal, non-bipartite, even delta"],
    [5, 11, 12,   1, 4, "Case 2a - Antipodal, non-bipartite, odd delta"], # Theorem 7.8.3

    [8, 17, 18,   1, 7, "Case 2a - Antipodal, non-bipartite, even delta"],
    # Note 100 is used as infinity, since 100 > 3 * delta in these cases.
    [3,  7, 10, 100, 0, "Case 1 - bipartite, non-antipodal, odd delta"], # Example 3, Corollary 6.10
    [4,  9, 12, 100, 0, 'Case 1 - bipartite, non-antipodal, even delta'],
    [6, 13, 16, 100, 0, 'Case 1 - bipartite, non-antipodal, even delta'],
    [4,  9, 10, 100, 0, 'Case 1 - bipartite, antipodal, even delta'], # Figure 12

    [5, 11, 12, 100, 0, 'Case 1 - bipartite, antipodal, odd delta'], # Theorem 7.8.2
    [6, 13, 14, 100, 0, 'Case 1 - bipartite, antipodal, even delta'], # Theorem 7.8.4
    [7, 15, 16, 100, 0, 'Case 1 - bipartite, antipodal, odd delta'],
    [8, 17, 18, 100, 0, 'Case 1 - bipartite, antipodal, even delta']
    ]

def graph_sample(n, parameters):
    """Return a premade graph (usually) of diameter delta.

    n should be chosen in range(1,15)
    The graphs [1,5,6,7,12,13] usually can't be completed, but the
    algorithm will try its best.
    """
    delta = parameters[0]
    number_of_vertices = delta + 1
    M = min_magic_number(parameters)

    G_sample_0 = Graph(1)

    # This one is a cycle with increasing distances and a disjoint point.
    # Don't expect this to be completed!
    G_sample_1 = Graph()
    for i in range(0,number_of_vertices):
        G_sample_1.add_edge((i,(i+1)%number_of_vertices, str(1 + (i%(number_of_vertices-1)))))
    G_sample_1.add_vertex(number_of_vertices+1)

    # This is a cycle of all 1s with delta many vertices.
    G_sample_2 = Graph()
    for i in range(1,number_of_vertices):
        G_sample_2.add_edge((i,(i+1)%number_of_vertices, "1"))

    # This is a path of 1s with delta+1 vertices.
    G_sample_3 = Graph()
    for i in range(number_of_vertices+1):
        G_sample_3.add_edge((i,(i+1), "1"))

    # This is a disjoint union of edges labelled 1
    G_sample_4 = Graph()
    for i in range(number_of_vertices // 2):
        G_sample_4.add_edge((2*i, 2*i + 1, "1"))

    # This a 1115 cycle from Figure 6
    # Don't expect this to be completed!
    G_sample_5 = Graph()
    for i in range(3):
        G_sample_5.add_edge((i, i + 1, "1"))
    G_sample_5.add_edge((3, 0, "5"))

    # This a 11555 cycle from Figure 7
    # Don't expect this to be completed!
    G_sample_6 = Graph()
    for i in range(2):
        G_sample_6.add_edge((i, i + 1, "1"))
    for i in range(2,5):
        G_sample_6.add_edge((i, (i + 1)%5, "5"))

    # This a 122 triangle from Figure 10
    # Don't expect this to be completed!
    G_sample_7 = Graph()
    G_sample_7.add_edges([(0, 1, "1"), (1, 2, "2"), (2, 0, "2")])

    # This is a disjoint union of edges labelled delta = 3 from Example 3
    G_sample_8 = Graph()
    G_sample_8.add_edges([(0, 1, "3"), (2, 3, "3")])

    # This is a (not disjoint) union of antipodal quadralaterals from Figure 12
    G_sample_9 = Graph()
    G_sample_9.add_edges([(1, 4, "1"), (2, 3, "1")])
    G_sample_9.add_edges([(3, 6, "2"), (4, 5, "2"), (3, 5, "2"), (4, 6, "2")])
    G_sample_9.add_edges([(1, 3, "3"), (2, 4, "3")])
    G_sample_9.add_edges([(1, 2, "4"), (3, 4, "4"), (5, 6, "4")])

    # From Corollary 6.10
    G_sample_10 = Graph()
    G_sample_10.add_vertex(2)
    G_sample_10.add_edge((0, 1, "1"))

    # From Theorem 7.8.2
    G_sample_11 = Graph()
    G_sample_11.add_vertex(2)
    G_sample_11.add_edge((0, 1, str(delta)))

    # From Theorem 7.8.3
    # Don't expect this to be completed!
    G_sample_12 = Graph()
    G_sample_12.add_edge((0, 1, str(delta)))
    G_sample_12.add_edges([(0, 2, str(M)), (1, 3, str(M)), (2, 3, str(M)), (2, 4, str(M)), (3, 4, str(M))])
    G_sample_12.add_edges([(0, 3, str(M+1)), (1, 2, str(M+1))])

    # From Theorem 7.8.4
    # Don't expect this to be completed!
    G_sample_13 = Graph()
    G_sample_13.add_edge((0, 1, str(delta)))
    G_sample_13.add_edges([(1, 2, str(delta // 2)), (0, 2, str(delta // 2))])
    G_sample_13.add_edge((2, 3, '1'))

    # A single edge
    G_sample_14 = Graph()
    G_sample_14.add_edge((0,1, "1"))

    graph_samples = [G_sample_0, G_sample_1, G_sample_2, G_sample_3,
                     G_sample_4, G_sample_5, G_sample_6, G_sample_7,
                     G_sample_8, G_sample_9, G_sample_10, G_sample_11,
                     G_sample_12, G_sample_13, G_sample_14]

    return graph_samples[n%len(graph_samples)]

#######################################
### 8.1 CHOOSE YOUR PARAMETERS HERE ################################################
#######################################

parameters_start = parameters_sample[15]

###########################################
### 8.2 CHOOSE YOUR STARTING GRAPH HERE ###########################################
###########################################

sample_graph_number = 1

G_start = graph_sample(sample_graph_number, parameters_start)

############################################################
### 8.3 Do you want all steps displayed? (True or False) ##########################
############################################################

display_all_steps_start = True
display_all_recursive_steps_start = True

###############################
### 8.4 The function call ##########################
###############################

display_graph_completion(G_start, parameters_start, display_all_steps_start, display_all_recursive_steps_start)
