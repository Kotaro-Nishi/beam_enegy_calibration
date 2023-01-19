import graphviz

g = graphviz.Digraph(format='png', filename='NonGUI_function')

#ModularEvaluation_NonGUI
with g.subgraph(name="cluster_ModularEvaluation_NonGUI") as c:
    c.attr(color='darkgreen')
    c.node("M1","per_iteation")
    c.node("M2","diffraction_dataset")
    c.node("M3","objective_diffraction")
    c.node("M4","MODFIT")

#UnDi
with g.subgraph(name="cluster_UnDi") as c:
    c.attr(color='red')
    c.node("U1","CoditionAmpl")
    c.node("U2","Conditionphase")
    c.node("U3","amplitudeonaperture_X")
    c.node("U4","phaseonaperture_X")
    c.node("U5","amplitudeonaperture")
    c.node("U6","phaseonaperture")
    c.node("U7","plotaperture")
    c.node("U8","superamplitudeInitialFit")
    c.node("U9","superamplitudeInitialFit_X")
    c.node("U10","plainintensity")
    c.node("U11","diffractionintensityInitialFit_X")
    c.node("U12","diffractionintensityInitialFit")

#AUndu
with g.subgraph(name="cluster_AUndu") as c:
    c.attr(color='blue')
    c.node("A1","Amlitude")
    c.node("A2","LamU")
    c.node("A3","Lu")

#IFunc
with g.subgraph(name="cluster_IFunc") as c:
    c.attr(color='black')
    c.node("I1","convolveplain")
    c.node("I2","convolvefresnel")
    c.node("I3","ElectronSourcefunction")
    c.node("I4","ElectronSourcefunctionSet")
    c.node("I5","squaredphase")

#plot_modular_resultsNonGUI
g.attr(color='yellow')
g.node("P","plot_modular_results_NonGUI")


####### 
#EDGES#
#######
g.edge("P","M4")
g.edge("U11","P")
g.edge("U10","P")
g.edge("U9","P")

g.edge("M3","M4")
g.edge("M1","M4")
g.edge("M2","M3")
g.edge("U11","M3")

g.edge("A1","U3")
g.edge("U1","U3")
g.edge("A2","U3")
g.edge("A2","U4")
g.edge("A3","U4")
g.edge("A2","U5")
g.edge("A1","U5")
g.edge("A2","U6")
g.edge("A3","U6")
g.edge("U3","U9")
g.edge("U6","U9")
g.edge("U9","U10")
g.edge("U10","U11")  
g.edge("I3","U11") 
g.edge("I5","U11") 
g.edge("I2","U11") 
g.edge("U8","U12")
g.edge("I3","U12")
g.edge("I5","U12")
g.edge("I2","U12")


g.render(directory="nishi",view =True)


h = graphviz.Digraph(format='png', filename='GUI_function')

with h.subgraph(name="cluster_Eval_Modular.py") as c:
    c.attr(color='purple')
    c.node("E1","load_button_clicked")
    c.node("E2","plot_button_clicked")
    c.node("E3","fit_button_clicked")
    c.node("E4","selectfile")

with h.subgraph(name="cluster_Modular_functions.py") as c:
    c.attr(color='darkgreen')
    c.node("M1","per_iteation")
    c.node("M2","diffraction_dataset")
    c.node("M3","objective_diffraction")
    c.node("M4","MODFIT")

with h.subgraph(name="cluster_Modular_functions.py") as c:
    c.attr(color='darkgreen')
    c.node("M1","per_iteation")
    c.node("M2","diffraction_dataset")
    c.node("M3","objective_diffraction")
    c.node("M4","MODFIT")

with h.subgraph(name="cluster_UnDi") as c:
    c.attr(color='red')
    c.node("U1","CoditionAmpl")
    c.node("U2","Conditionphase")
    c.node("U3","amplitudeonaperture_X")
    c.node("U4","phaseonaperture_X")
    c.node("U5","amplitudeonaperture")
    c.node("U6","phaseonaperture")
    c.node("U7","plotaperture")
    c.node("U8","superamplitudeInitialFit")
    c.node("U9","superamplitudeInitialFit_X")
    c.node("U10","plainintensity")
    c.node("U11","diffractionintensityInitialFit_X")
    c.node("U12","diffractionintensityInitialFit")

with h.subgraph(name="cluster_IFunc") as c:
    c.attr(color='black')
    c.node("I1","convolveplain")
    c.node("I2","convolvefresnel")
    c.node("I3","ElectronSourcefunction")
    c.node("I4","ElectronSourcefunctionSet")
    c.node("I5","squaredphase")

with h.subgraph(name="cluster_AUndu") as c:
    c.attr(color='blue')
    c.node("A1","Amlitude")
    c.node("A2","LamU")
    c.node("A3","Lu")

with h.subgraph(name="cluster_Ana") as c:
    c.node("N1","Ana")

#EDGES
h.edge("N1","E1")



h.render(directory="nishi",view=True)