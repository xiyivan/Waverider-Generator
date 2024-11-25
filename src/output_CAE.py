from taylor_maccoll_sol import Taylor_Maccoll as tm
from TE_Formation import TEG
from streamline_tracing import TRACE

class CAE:
    def __init__(self) -> None:
        self.TEG = TEG()
        self.trace = TRACE()

    def output_bp(self, a, b, c, Rs, L, N, N_l, N_up, Vr_i, V_theta_i):

        X_b, Y_b, Z_b, X_p, Y_p, Z_p = self.trace.projection_module(a, b, c, Rs, L, N)
        m_disp_carte_x, m_disp_carte_y, m_disp_carte_z = self.trace.tracing_module(L, N, N_l, N_up, Vr_i, V_theta_i)

        X_p[0] = L
        X_p[-1] = L  # This is to eliminate the uncertainty arised from significant figures.
                     # The uncertainty is smaller than 1e-15, however for better external CAD modeling, 
                     # the start and end values are still replaced with L
        with open('Leading Edge.txt', 'w') as file:
            for i in range(len(X_p)):
                file.write(f"{X_p[i]}\t{Y_p[i]}\t{Z_p[i]}\n")

        X_b_out = []
        for i in range(0, N):
            X_b_out.append(X_b)

        with open('Trailing Edge.txt', 'w') as file:
            for i in range(len(Y_b)):
                file.write(f"{X_b_out[i]}\t{Y_b[i]}\t{Z_b[i]}\n")

        m_disp_carte_x_out = []
        for i in range(0, len(m_disp_carte_y)):
            m_disp_carte_x_out.append(L)

        with open('Lower Surface Base Plane Curve.txt', 'w') as file:
            for i in range(len(m_disp_carte_x)):
                file.write(f"{m_disp_carte_x_out[i]}\t{m_disp_carte_y[i]}\t{m_disp_carte_z[i]}\n")

    def output_lower_surface(self):
        for idx, curve_data in enumerate(self.trace.lower_surface):

            with open(f"{idx} Lower Surface Line.txt", "w") as file:
                for point in curve_data["curve"]:
                    file.write(f"{point[0]}\t{point[1]}\t{point[2]}\n")

            with open(f"{idx} Mirrored Lower Surface Line.txt", "w") as file:
                for point in curve_data["mirrored_curve"]:
                    file.write(f"{point[0]}\t{point[1]}\t{point[2]}\n")
