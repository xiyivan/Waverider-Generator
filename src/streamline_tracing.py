import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev

from taylor_maccoll_sol import Taylor_Maccoll as tm
from TE_Formation import TEG

class TRACE:
    def __init__(self) -> None:
        # Attribute
        self.lower_surface = []

        # Class Object
        self.tm_i = tm()
        self.TEG = TEG()

    def projection_module(self, a, b, c, Rs, L, N):

        curve = TEG.te_plot(self, a, b, c, Rs, L, N)
        y_plot = curve[0]
        z_plot = curve[1]
        self.x_plot_break = curve[2]
        self.y_plot_break = curve[3]
        self.z_plot_break = curve[4]

        #   Leading Edge
        self.X_p = L * (1-(Rs-np.sqrt(y_plot**2 + z_plot**2))/Rs)
        self.Y_p = y_plot
        self.Z_p = z_plot

        #   Trailing Edge
        self.X_b = L
        self.Y_b = y_plot
        self.Z_b = z_plot

        return self.X_b, self.Y_b, self.Z_b, self.X_p, self.Y_p, self.Z_p

    def tracing_module(self, L, N, N_l, N_up, Vr_i, V_theta_i):

        # Initialize Data Storage Matrices
        baseplane_carte_x = []
        baseplane_carte_y = []
        baseplane_carte_z = []
        baseplane_carte_z_mir = []
        baseplane_carte_y_mir = []

        ax = plt.axes(projection = '3d')

        for temp in range(0, N_l):

            r_i = np.sqrt(self.x_plot_break**2 + (self.z_plot_break)**2 + self.y_plot_break**2)

            idx = [0]

            for i in range(1, N_l+1):
                j = len(r_i) % N_l
                if j == 0:
                    idx.append(int(i*(len(r_i)/N_l)-1))
                else:
                    idx.append(int(i*((len(r_i)-j)/N_l)-1))

            #print(len(r_i))
            #print(idx)

            phi = -np.atan(self.y_plot_break[idx[temp]]/self.z_plot_break[idx[temp]])
            print(np.degrees(phi))

            def event_cr2(theta, S):
                return S[1]
            event_cr2.terminal = True

            r_i_s = np.sqrt(self.x_plot_break**2 + (self.z_plot_break)**2 + self.y_plot_break**2)

            theta_i_s = np.atan(np.sqrt((self.z_plot_break)**2 + self.y_plot_break**2)/self.x_plot_break)

            #   Conical Shock Solution, Loop Iteration, Data Stored for Euler Integration
            thetas = [np.abs(theta_i_s[idx[temp]]), 1e-08]
            theta_range = np.linspace(np.abs(theta_i_s[idx[temp]]), 1e-08, 1000)
            sol2 = self.tm_i.tracing_solver(Vr_i, V_theta_i, thetas, theta_range)

            #   re-dimensionize V'r, V'theta not needed, note V'r/V'theta = Vr/Vtheta = S2; define stream eqn
            S2 = []
            for i in range(len(sol2.t)):
                S2.append(sol2.y[0][i]/sol2.y[1][i])

            Vr_Vtheta_ratio = np.array(S2)
            thet_array = np.array(sol2.t)

            print(np.degrees(thet_array))

            #   Uncomment to Check Case: len(Vr_Vtheta_ratio) should equal len(thet_array)
            '''
            print(Vr_Vtheta_ratio)
            print(len(Vr_Vtheta_ratio))
            print(thet_array)
            print(len(thet_array))
            '''

            #   Euler Integration With Different Initial Condition, Same Stop Criterion, Loop Iteration
            r_march = [r_i_s[idx[temp]]]
            thet_march = [thet_array[0]]

            for i in range(1, len(thet_array)):
                d_theta = thet_array[i] - thet_array[i-1]
                r_i_s = r_march[-1]
                Vr_Vtheta = Vr_Vtheta_ratio[i-1]

                dr_dtheta = r_i_s * Vr_Vtheta
                r_n = r_i_s + (dr_dtheta * d_theta)

                # Last Trace Point Interpolation, Convert to cartesian, Transform, Then Convert Back

                if r_i_s*np.cos(thet_array[i-1]) - L < 0 <= r_n*np.cos(thet_array[i]) - L:

                    z_itp = np.cos(phi)*((r_n*np.sin(thet_array[i])-r_i_s*np.sin(thet_array[i-1]))/(r_n*np.cos(thet_array[i])-r_i_s*np.cos(thet_array[i-1])))*L + \
                            np.cos(phi)*(r_i_s*np.sin(thet_array[i-1])-r_i_s*np.cos(thet_array[i-1])*(r_n*np.sin(thet_array[i])-r_i_s*np.sin(thet_array[i-1]))/(r_n*np.cos(thet_array[i])-r_i_s*np.cos(thet_array[i-1])))

                    r_break = np.sqrt((z_itp/np.cos(phi))**2 + L**2)
                    theta_break = np.atan((z_itp/np.cos(phi))/L)

                    r_march.append(r_break)
                    thet_march.append(theta_break)

                    print(f"Euler Integration Stop at r_march[{i}] = {r_n}, with interpolated theta = {np.degrees(theta_break)}")
                    break
                else:
                    r_march.append(r_n)
                    thet_march.append(thet_array[i])

            carte_x = []
            carte_z = []
            carte_y = []
            carte_y_mir = []

            # Store Data for Plotting Traced Streamline and Output to CAD Software
            # Translate Streamlines to Inclined Osculating Planes

            for i in range(len(r_march)):
                theta_m_deg = np.degrees(thet_march[i])
                r_m = r_march[i]
                print(f"Theta: {theta_m_deg: .4f} degrees, Radius: {r_m: .4f}")

                carte_z.append(-1*r_m*np.sin(thet_march[i])*np.cos(phi))

                carte_x.append(r_m*np.cos(thet_march[i]))
            
                carte_y.append(r_m*np.sin(thet_march[i])*np.sin(phi))
                carte_y_mir.append(-r_m*np.sin(thet_march[i])*np.sin(phi))

                print(f"Z: {carte_z[i]: .4f}, X: {carte_x[i]: .4f}, Y: {carte_y[i]: .4f}, Phi: {np.degrees(phi)}")

            baseplane_carte_z.append(carte_z[-1])
            baseplane_carte_z_mir.append(carte_z[-1])

            baseplane_carte_y.append(carte_y[-1])
            baseplane_carte_y_mir.append(-carte_y[-1])

            baseplane_carte_x.append(carte_x[-1])

            ax.plot(carte_x, carte_y, carte_z, color = 'b')
            ax.plot(carte_x, carte_y_mir, carte_z, color = 'b')

            carte_x[-1] = L
            crv = []
            crv_mir = []

            for i in range(len(carte_x)):

                crv_cord = (carte_x[i], carte_y[i], carte_z[i])
                crv.append(crv_cord)

                # Append mirrored points (x, -y, z) for the mirrored curve
                crv_cord_mir = (carte_x[i], carte_y_mir[i], carte_z[i])
                crv_mir.append(crv_cord_mir)

            ls_curve_dat = {
                "curve": crv,
                "mirrored_curve": crv_mir
            }
            self.lower_surface.append(ls_curve_dat)

        baseplane_carte_z.append(self.z_plot_break[-1])
        baseplane_carte_z_mir.append(self.z_plot_break[-1])
        baseplane_carte_z.reverse()
        baseplane_carte_z_mir.pop(0)
        m_disp_carte_z = baseplane_carte_z + baseplane_carte_z_mir

        baseplane_carte_y.append(self.y_plot_break[-1])
        baseplane_carte_y_mir.append(-self.y_plot_break[-1])
        baseplane_carte_y.reverse()
        baseplane_carte_y_mir.pop(0)
        m_disp_carte_y = baseplane_carte_y + baseplane_carte_y_mir

        m_disp_carte_x = baseplane_carte_x + baseplane_carte_x
        m_disp_carte_x.append(self.x_plot_break[-1])

        #print(f'Length of X is {len(m_disp_carte_x)}, Y is {len(m_disp_carte_y)}, Z is {len(m_disp_carte_z)}')
        tck, u = splprep([m_disp_carte_x, m_disp_carte_y, m_disp_carte_z], s=0)
        u_fine = np.linspace(0, 1, 1000)
        x_smooth, y_smooth, z_smooth = splev(u_fine, tck)

        # Plot the Waverider Geometry
        ax.plot(x_smooth, y_smooth, z_smooth, color='r', lw=2, label='Interpolated Lower Surface Curve')
        ax.plot(self.X_p, self.Y_p, self.Z_p, color = 'k', label='Leading Edge')
        ax.plot(self.X_b, self.Y_b, self.Z_b, color = 'y', label='Trailing Edge')

        # Make the Two Lines Below to Comments if Only Want to Display Waverider Geometry
        #ax.plot_surface(X, Y, Z, alpha = 0.5, cmap = 'Greys', edgecolor = 'none', label='Conical Shock Wave')
        #ax.plot_surface(X_c, Y_c, Z_c, alpha = 0.1, label='Base Cone')

        for temp2 in range(0, N_up):

            idx_up = [0]

            for i in range(1, N_up):
                j = len(r_i) % N_up
                if j == 0:
                    idx_up.append(int(i*(len(r_i)/N_up)-1))
                else:
                    idx_up.append(int(i*((len(r_i)-j)/N_up)-1))

            x_u = np.linspace(self.x_plot_break[idx_up[temp2]], L, N)
            y_u = self.y_plot_break[idx_up[temp2]]
            z_u = self.z_plot_break[idx_up[temp2]]

            plt.plot(x_u, y_u, z_u, color = 'c')
            plt.plot(x_u, -y_u, z_u, color = 'c')

        plt.legend(loc='best')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.axis('equal')

        return m_disp_carte_x, m_disp_carte_y, m_disp_carte_z