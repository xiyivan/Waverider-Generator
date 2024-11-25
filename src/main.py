import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from scipy.interpolate import splprep, splev

from oblique_shock import Oblique_Shock
import user_input as UI
from taylor_maccoll_sol import Taylor_Maccoll
from TE_Formation import TEG
from streamline_tracing import TRACE
from output_CAE import CAE

class Main:
    def __init__(self) -> None:
       
        # Attribute
        self.M1 = UI.M1
        self.gamma = UI.gamma
        self.beta = UI.beta
        self.beta_rad = UI.beta_rad
        self.a = UI.a
        self.b = UI.b
        self.c = UI.c
        self.Rs = UI.Rs
        self.L = UI.L
        self.N = UI.N
        self.N_l = UI.N_l
        self.N_up = UI.N_up

        # Class Object
        self.os = Oblique_Shock()
        self.tm = Taylor_Maccoll()
        self.TEG = TEG()
        self.trace = TRACE()
        self.cae = CAE()

    def run(self):
        #    Using Oblique Shock Relation Obtain Initial Conditions
        self.os.sub_1(self.M1, self.gamma, self.beta)
        BC = self.os.initial_nondimensioned_conditions(self.M1, self.gamma, self.beta)
        Vr_i = BC[0]
        V_theta_i = BC[1]

        #    Solving Non-Dimensioned Taylor-Maccoll Equations
        self.tm.solver1_result_visualize(self.beta_rad, Vr_i, V_theta_i)

        #    TE, intake 3 const, specify your UDF in UserInput Under src/
        #    The second line plots TE and conical shock on baseplane
        self.TEG.te_plot(self.a, self.b, self.c, self.Rs, self.L, self.N)
        self.TEG.baseplane_visualize(self.L, self.Rs, self.N, self.beta_rad, Vr_i, V_theta_i)
        plt.show()

        #    LE Forward Projection onto Shock Cone
        self.trace.projection_module(self.a, self.b, self.c, self.Rs, self.L, self.N)

        #    Streamline Tracing Marching from Shock Surface to Baseplane
        self.trace.tracing_module(self.L, self.N, self.N_l, self.N_up, Vr_i, V_theta_i)
        plt.show()

        '''
            | OUTPUT MODULE

            | Format: .txt file of Key Curve Coordinates for Curve Import in CAD/CAE Software
            | Unit:   MKS     [Meters,             Kilograms,            Seconds]
            | Comment out the following line if you don't want any outputs to your device yet

           \  /  CAUTION: Those line will write N_l+N_up+2 number of .txt file to your device
            \/
        ''' 
        #self.cae.output_bp(self.a, self.b, self.c, self.Rs, self.L, self.N, self.N_l, self.N_up, Vr_i, V_theta_i)
        #self.cae.output_lower_surface()
    
main = Main()
main.run()
