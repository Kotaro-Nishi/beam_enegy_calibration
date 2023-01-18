import graphviz

g = graphviz.Digraph(format='png', filename='GUI')

g.edge("MoFu","Eval_Modular")
g.edge("AUndu","IFunc")
g.edge("UnDi","ModularEvaluation")
g.edge("IFunc","ModularEvaluation")
g.edge("AUndu","ModularEvaluation")
g.edge("UnDi","MoFu")
g.edge("IFunc","MoFu")
g.edge("AUndu","MoFu")
g.edge("UnDi","plot_modular_results")
g.edge("IFunc","UnDi")
g.edge("AUndu","UnDi")
g.edge("tkFIT","Eval_Modular")

g.render(directory="nishi",view =True)

h = graphviz.Digraph(format='png',filename ='NON_GUI')

h.edge("AUndu","IFunc")
h.edge("UnDi","ModularEvaluation_NonGUI")
h.edge("IFunc","ModularEvaluation_NonGUI")
h.edge("AUndu","ModularEvaluation_NonGUI")
h.edge("AUndu","plot")
h.edge("UnDi","plot")
h.edge("plot","ModularEvaluation_NonGUI")

h.render(directory="nishi",view =True)
