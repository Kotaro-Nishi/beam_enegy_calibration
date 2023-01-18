import graphviz

g = graphviz.Digraph(format='png', filename='GUI')

g.node("A","Alice")
g.node("B","Bob")
g.edges(["AB"])

g.render(directory="nishi",view =True)