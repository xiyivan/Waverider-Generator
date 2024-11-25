import numpy as np
from scipy.integrate import solve_ivp
from user_input import gamma

class Taylor_Maccoll:
    def __init__(self) -> None:
        pass

    def TM_eqn(self, theta, S):

        A, A_pr = S  # A = Vr, A_pr = dVr/dtheta

        tep1 = (A * A_pr**2 - (gamma - 1) / 2 * (1 - A**2 - A_pr**2) * (2 * A + A_pr / np.tan(theta)))
        tep2 = ((gamma - 1) / 2 * (1 - A**2 - A_pr**2) - A_pr**2)
        A_2nd = tep1 / tep2

        func = [
            A_pr,          # A_pr = dVr/dtheta
            A_2nd,         # A_2nd = d2Vr/dtheta2 = J.D Anderson Modern Compressible Flow Eqn 10.15
        ]
        return func

    def solver(self, beta_rad, Vr_i, V_theta_i):

        def event_cr(theta, S):
            return S[1]
        event_cr.terminal = True

        thetas = [beta_rad, 1e-08]
        A_0 = Vr_i
        A_pr0 = V_theta_i
        S_0 = (A_0, A_pr0)
        sol = solve_ivp(self.TM_eqn, thetas, y0=S_0, method='RK45', events=event_cr, rtol=1e-08, atol=1e-10)

        return sol

    def solver1_result_visualize(self, beta_rad, Vr_i, V_theta_i):

        sol = self.solver(beta_rad, Vr_i, V_theta_i)

        for i in range(len(sol.t)):
            theta_deg = np.degrees(sol.t[i])
            print(f"Theta: {theta_deg:.2f} degrees, Vr': {sol.y[0][i]:.4f}, Vtheta': {sol.y[1][i]:.4f}")

        cone = sol.t[-1]
        print(f'The Base Cone is {np.degrees(cone)} degrees')

        return cone
        
    def tracing_solver(self, Vr_i, V_theta_i, thetas, theta_range):

        def event_cr2(theta, S):
            return S[1]
        event_cr2.terminal = True

        A_02 = Vr_i
        A_pr02 = V_theta_i
        S_02 = (A_02, A_pr02)
        sol2 = solve_ivp(self.TM_eqn, thetas, y0=S_02, method='RK45', events=event_cr2, t_eval=theta_range, rtol=1e-08, atol=1e-10)

        return sol2

