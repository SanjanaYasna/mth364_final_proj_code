#get data needed for bifurcation of $NF_1 S^{B+}_{A+}$ 
#in terms of B* from different values of I

from gekko import GEKKO
import numpy as np
from torch.func import jacrev
from torch import tensor

#csv sholuld have:
#I, A, B, trace, determinant, string_of_classification 

with open("NF1_I_irreversible.csv", "w") as f:
    f.write("I,A,B,trace,determinant,string_of_classification,stable\n")
f.close() 


values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
m = GEKKO(remote =False)   

m.options.SOLVER=3
m.options.MAX_ITER=50000
m.options.COLDSTART=1 

#function for NF1, for computing jacobian with
def f_w_params(A, B,
               I, k_I, k_A, k_B, k_1, k_2,
               k_mI, k_mA, k_mB, k_m1, k_m2,
               delta_1, delta_2):
    return (
        ((k_I * I * (1 - A)) / (k_mI + 1 - A) ) +( (k_A * A * (1- A)) / (k_mA + 1 - A) ) + ((-k_2*B*A) / (k_m2 + A)) - (delta_1 * A) ,
        (( k_B * B * (1- B)) / (k_mB+ 1 - B))+((k_1*A*(1-B)) / (k_m1 + 1 - B) )- (delta_2 * B)
    )
    
#iterate between values of I from 0 to 1 with step 0.01
I_values =np.arange(0, 1.0, 0.01)

for I in I_values:
    #unique_solutions tracks solutions for a given value of I
    unique_solutions = [] 
    #gekko does grid search over solutions (9 x 9 ics)
    for i in values:
        for j in values:
            m = GEKKO(remote =False)   
            #all params below are the same except I
            I = I
            #convert to gekko param
            I = m.Param(value=I)
            
            #all other params that will be held constant
            #redundant, but need to convert to gekko params for each iteration...doesn't work otherwise 
            k_I = 0.1867 
            k_mI = 0.7801
            k_A = 0.9833 
            k_mA = 0.4083
            k_B = 0.8252
            k_mB = 0.5157
            k_1 = 0.3132
            k_m1= 0.5770 
            #k_2 =0.3700 #REVERSIBLE
            k_2 = 0.2636 # is IRREVERSIBLE
            k_m2 = 0.1090 
            delta_1 = 0.0697 
            delta_2 = 0.0226

            #convert them to gekko-type paramters (gekko is solver)
            k_I = m.Param(value=k_I)
            k_mI = m.Param(value=k_mI)
            k_A = m.Param(value=k_A)
            k_mA = m.Param(value=k_mA)
            k_B = m.Param(value=k_B)
            k_mB = m.Param(value=k_mB)
            k_1 = m.Param(value=k_1)
            k_m1 = m.Param(value=k_m1)
            k_2 = m.Param(value=k_2)
            k_m2 = m.Param(value=k_m2)
            delta_1 = m.Param(value=delta_1)
            delta_2 = m.Param(value=delta_2)
            

            #params for solver
            m.options.SOLVER=3
            m.options.MAX_ITER=100
            m.options.COLDSTART=1 
            m.options.OTOl =    1e-20
            
            #gekko variables for steady state
            A = m.Var(value=i)
            A.lower = 0
            A.upper = 1
            B = m.Var(value=j)
            B.lower = 0
            B.upper = 1
            
            m.Equations([
                #dA/dt is a doozy
                ((k_I * I * (1 - A)) / (k_mI + 1 - A) ) +( (k_A * A * (1- A)) / (k_mA + 1 - A) ) + ((-k_2*B*A) / (k_m2 + A)) - (delta_1 * A) == 0,
                #dB/dt
                (( k_B * B * (1- B)) / (k_mB+ 1 - B))+((k_1*A*(1-B)) / (k_m1 + 1 - B) )- (delta_2 * B) == 0 ])
            
            try:
                m.solve(disp=False)
                #round to 3 numbers
                A.value[0] = round(A.value[0], 4)
                B.value[0] = round(B.value[0], 4)
                #add solution only if it hasn't already been recorded 
                if (A.value[0], B.value[0]) not in unique_solutions:
                    unique_solutions.append((A.value[0], B.value[0]))
            except: 
                #exception happens when solution fails, but that's expected
                pass
            
    #convert variables to tensors:
    I = tensor(I)
    k_I = tensor(k_I)
    k_A = tensor(k_A)
    k_B = tensor(k_B)
    k_1 = tensor(k_1)
    k_2 = tensor(k_2)
    k_mI = tensor(k_mI)
    k_mA = tensor(k_mA)
    k_mB = tensor(k_mB)
    k_m1 = tensor(k_m1)
    k_m2 = tensor(k_m2)
    delta_1 = tensor(delta_1)
    delta_2 = tensor(delta_2)
    #now that we have the solutions for the vlaue of I, compute stabilities of the solutions
    #iterate over solution
    
    for solution in unique_solutions:
        A = tensor(solution[0])
        B = tensor(solution[1])
        #jacobian with argnums specifying which paramters are arguments to take jacobianwith respect to
        J = jacrev(f_w_params, argnums=(0, 1))(A, B, I, k_I, k_A, k_B, k_1, k_2,
            k_mI, k_mA, k_mB, k_m1, k_m2,
            delta_1, delta_2)
        J = np.array(J)
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        #if determinant is nggative 
        #compute trace
        trace = np.trace(J)
        classification = ""
        stable = False
        #if trace squared is less than 4 times the determinant
        if trace**2 < 4*det:
            classification = "stable spiral node sink"
            stable = True
        elif det > 0 and trace < 0:
            classification = "stable node sink" 
            stable = True
        elif det < 0:
            classification = "unstable saddle node"
            stable = False
        elif trace**2 > 4*det and trace > 0:
            classification = "unstable node source"
            stable = False
        else:
            classification = "unstable spiral source"
            stable = False
        #write to csv
        with open("NF1_I_irreversible.csv", "a") as f:
            #cnovert I to number
            I = float(I)
            #get trace
            trace = trace[0]
            det = det[0]    
            f.write(f"{I},{A},{B},{trace},{det},{classification},{stable}\n")
        f.close()
        print(f"I: {I}, A: {A}, B: {B}, trace: {trace}, determinant: {det}, classification: {classification}, stable: {stable}")