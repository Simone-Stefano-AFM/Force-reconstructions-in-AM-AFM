# Matrix FRM for AM-AFM
This is a Python code to reconstruct tip-sample forces in Amplitude Modulation Atomic Force Microscopy.

---

## Description and How to Use
The two algorithms are interconnected. Please run Holscher_FRM.py before Matrix.py, as the former generates k_Holscher.txt, which contains the necessary information for the second algorithm to extract the force. 
k_Holscher.txt is a two-column dataset where the first column represents the minimum tip-sample distance, and the second column contains basically the virial of the tip-sample force.
Holscher_FRM.py allows to reconstruct the tip-sample force with the original formula developed by Holscher (see doi.org/10.1063/1.2355437).
Matrix.py, instead, is a mathematical method that has been developed to overcome limitations of the Holscher algorithm.

---

## Additional Information and How to Cite
Additional information on the tool (Matrix.py) can be found in the original article "Quantification of solvation forces with amplitude modulation AFM", at the link https://www.sciencedirect.com/science/article/pii/S0021979725001456#da005.
If you use the matrix-based algorithm, please cite the article above (doi.org/10.1016/j.jcis.2025.01.131).
