import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from taylor_maccoll_sol import Taylor_Maccoll as tm
import user_input as UI

class TEG:
    def __init__(self) -> None:
        # Class Object
        self.tm_i = tm()

    def te_plot(self, a, b, c, Rs, L, N):

        def equation(y):
            return np.sqrt(y**2 + UI.z(y)**2) - Rs

        y_up = fsolve(equation, 2.0)[0]   # Solve for the Upper Bound (Coordinate Which Satisfies sqrt(z**2 + y**2) = Rs)

        self.y_plot = np.linspace(-y_up, y_up, N)
        self.z_plot = UI.z(self.y_plot)

        #   This is used as later Euler Integration initial value, and check case to break loop
        y_plot_break = np.linspace(0, y_up, N)
        z_plot_break = UI.z(y_plot_break)
        x_plot_break = L * (1-(Rs-np.sqrt(y_plot_break**2 + z_plot_break**2))/Rs)

        #   Base Plane Conical Shock

        return self.y_plot, self.z_plot, x_plot_break, y_plot_break, z_plot_break
    
    def baseplane_visualize(self, L, Rs, N, beta_rad, Vr_i, V_theta_i):

        self.cone = self.tm_i.solver1_result_visualize(beta_rad, Vr_i, V_theta_i)

        y_conical = np.linspace(-Rs, Rs, N)
        z_conical = np.sqrt(Rs**2 - y_conical**2)
        z_conical_sym = -np.sqrt(Rs**2 - y_conical**2)

        r = np.linspace(0, Rs, N)
        theta_plot = np.linspace(0, 2*np.pi, N)
        R, THETA = np.meshgrid(r, theta_plot)
        X = R / np.tan(beta_rad)
        Y = R * np.cos(THETA)
        Z = R * np.sin(THETA)

        #   Base Cone Visualization
        R_c = L * np.tan(self.cone)

        Y_c = np.linspace(-R_c, R_c, N)
        Z_c = np.sqrt(R_c**2 - Y_c**2)
        Z_c_sym = -np.sqrt(R_c**2 - Y_c**2)

        #   Plot the Base Plane Conical Shock
        plt.plot(y_conical, z_conical, label=r'Conical Shock', color = 'r')
        plt.plot(y_conical, z_conical_sym, color = 'r')
        plt.plot(Y_c, Z_c, label=r'Base Cone', color = 'k')
        plt.plot(Y_c, Z_c_sym, color = 'k')

        #   Plot the Trailing Edge Function
        plt.legend(loc='best')
        plt.plot(self.y_plot, self.z_plot, label=r'Trailing Edge')
        plt.xlabel('y')
        plt.ylabel('z')
        plt.title(f'Base Plane')
        plt.axis('equal')
        plt.axis('equal')
        plt.grid(True)
        

    
