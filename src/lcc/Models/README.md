# Model Examples

## Ingalls, B. P. (2013). Mathematical modeling in systems biology: an introduction. Cambridge, Massachusetts: MIT Press.

### 1. Numerical Simulation (p.35)
<img src="Ingalls2013_Model2.18_NumericalSimulation_Model.png" height="50"> 
<img src="Ingalls2013_Model2.18_NumericalSimulation_Eqn.png" height="75">

    Initial Concentrations:     A = 0, B = 10
    Kinetic Constants:          k_1 = 9, k_-1 = 12, k_2 = 2
    Time interval:              t = 1/100

- [NumericalSimulation.py](Ingalls2013_Model2.18_NumericalSimulation.py)
- [NumericalSimulation.lpp](Ingalls2013_Model2.18_NumericalSimulation.lpp)


### 2. Full vs. Reduced (Michaelis Menten) Model (p.49)
<img src="Ingalls2013_Model3.2_MichaelisMenten_Model.png" height="50"> 
<img src="Ingalls2013_Model3.2_MichaelisMenten_Eqn1.png" height="150"> 
<img src="Ingalls2013_Model3.2_MichaelisMenten_Eqn2.png" height="50"> 
    
    Initial Concentrations:     S = 5, E = 1, C = 0, P = 0
    Kinetic Constants:          Full model      | k_1 = 30, k_-1 = 1, k_2 = 10
                                Reduced model   | KM = (krev1 + k2) / k1
    Time interval:              t = 1/500

- [MichaelisMentenKinetics.py](Ingalls2013_Model3.2_MichaelisMenten.py)
- [MichaelisMentenKinetics.lpp](Ingalls2013_Model3.2_MichaelisMenten.lpp)
- Examples:


### 3. Product Inhibition (p.78)
<img src="Ingalls2013_Model4.1_ProductInhibition_Model.png" height="100"> 
<img src="Ingalls2013_Model4.1_ProductInhibition_Eqn.png" height="100">

    Initial Concentrations:     S1 = 0, S2 = 0
    Kinetic Constants:          k1 = 20, k2 = 5, k3 = 5, k4 = 5, k5 = 2, K = 1, n = 4 

- [ProductInhibition.py](Ingalls2013_Model4.1_ProductInhibition.py)
- [ProductInhibition.lpp](Ingalls2013_Model4.1_ProductInhibition.lpp)
- Examples: 

### 4. Stability (p.82)
<img src="Ingalls2013_Model4.2_Stability_Model.png" height="100"> 
<img src="Ingalls2013_Model4.2_Stability_Eqn.png" height="100">

    Initial Concentrations:     S1 = 1, S2 = 3
    Kinetic Constants:          k1 = 20, k2 = 20, k3 = 5, k4 = 5, K1 = 1, K2 = 1, n1 = 4, n2 = 1 

- [ProductInhibition.py](Ingalls2013_Model4.2_Stability.py)
- [ProductInhibition.lpp](Ingalls2013_Model4.2_Stability.lpp)
- Examples: 







[comment]: <> (### Competitive Inhibition)

[comment]: <> (<img src="Ingalls2013_Model3.13x_CompetitiveInhibition_Model.png" height="100"> )

[comment]: <> (<img src="Ingalls2013_Model3.13x_CompetitiveInhibition_Eqn.png" height="50"> )

[comment]: <> (    Initial Concentrations:     S = [0, 1, ..., 100], E = 1, I = [0, 5, 10, 15])

[comment]: <> (    Kinetic Constants:          k_1 = 5, k_-1 = 1, k2 = 8, k3 = 2, k-3 = 1)

[comment]: <> (- [CompetitiveInhibition.py]&#40;Ingalls2013_Model3.13x_CompetitiveInhibition.py&#41;)

[comment]: <> ([comment]: <> &#40;- [AllostericRegulation.lpp]&#40;Ingalls2013_Model3.13x_CompetitiveInhibition.lpp&#41;&#41;)

[comment]: <> (- Examples: ibuprofen &#40;Nonsteroidal anti-inflammatory drug&#41;)


[comment]: <> (### Allosteric Regulation)

[comment]: <> (<img src="Ingalls2013_Model3.14_AllostericRegulation_Model.png" height="200"> )

[comment]: <> (<img src="Ingalls2013_Model3.14_AllostericRegulation_Eqn.png" height="50"> )

[comment]: <> (    Initial Concentrations:     S = [0, 1, ..., 50], E = 1, I = [0, 1.5, 3, 4.5])

[comment]: <> (    Kinetic Constants:          k_1 = 5, k_-1 = 1, k2 = 8, k3 = 2, k-3 = 1)

[comment]: <> (- [AllostericRegulation.py]&#40;Ingalls2013_Model3.14_AllostericRegulation.py&#41;)

[comment]: <> (- [AllostericRegulation.lpp]&#40;Ingalls2013_Model3.14_AllostericRegulation.lpp&#41;)

[comment]: <> (- Examples: benzodiazepines &#40;depressants&#41;)


[comment]: <> (### Cooperativity: Hill Function)

[comment]: <> (<img src="Ingalls2013_Model3.16_Cooperativity_Model.png" height="200"> )

[comment]: <> (<img src="Ingalls2013_Model3.16_Cooperativity_Eqn.png" height="50"> )

[comment]: <> (    Initial Concentrations:     X = [0, 1, ..., 200])

[comment]: <> (    Kinetic Constants:          K, n = [[5, 1], [20, 2], [45, 3], [80, 4]] )

[comment]: <> (- [Cooperativity.py]&#40;Ingalls2013_Model3.16_Cooperativity.py&#41;)

[comment]: <> ([comment]: <> &#40;- [Cooperativity.lpp]&#40;Ingalls2013_Model3.16_Cooperativity.lpp&#41;&#41;)

[comment]: <> (- Examples: Oxygen binding to Hemoglobin &#40;sigmoidal&#41; vs. Myoglobin &#40;hyperbolic&#41;)