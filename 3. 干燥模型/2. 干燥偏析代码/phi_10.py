import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def grap_phi_1(Pe1,Pe2,E,phi1_0,phi2_0,N = 10):
    # 参数定义
    # Pe1 = 0.18
    # Pe2 = 0.37
    # phi1_0 = np.linspace(0.101, 0.209, N)
    # phi2_0 = np.linspace(0.101, 0.209, N)
    
    
    if len(phi1_0) != len(phi2_0) and len(phi2_0) != N:
        raise ValueError("phi1,phi2,的初始形状必须跟时间空间矩阵形状一样。 就是要phi1，和phi2 的数据个数必须等于N,N默认为10")
        return None


    phi_max = 0.64
    # 空间和时间范围
    xi_min, xi_max = 0, 1
    tau_min, tau_max = 0.17, 1  # 避免 t 直达1
  
    num_xi = N
    num_tau = N
    # E = 1

    xi = np.linspace(xi_min, xi_max, num_xi)

    t = np.linspace(tau_min, tau_max, num_tau)



    def K(phi1, phi2):
        phi1_mean = np.mean(phi1)
        phi2_mean = np.mean(phi2)
        value = 1 - np.where(np.isfinite(phi1_mean) & np.isfinite(phi2_mean), phi1_mean + phi2_mean, 0)
        return value**6.55

    def Z(phi1, phi2):
        phi1_mean = np.mean(phi1)
        phi2_mean = np.mean(phi2)
        value = phi_max - phi1_mean - phi2_mean
        return phi_max / value
        # return 0.23


    def pde_system(t,y,Pe1,Pe2,xi):
        
        
        phi1 = y[:N]
        phi2 = y[N:]

        dphi1_dxi = np.gradient(phi1,xi)
        dphi2_dxi = np.gradient(phi2,xi)
        
        term1_phi1 = (K(phi1, phi2) * phi1 * ((phi1 * (1-phi1)) * (phi2/phi1) * (dphi1_dxi / dphi2_dxi) - phi1 * phi2 * (((Pe1 / Pe2)**3) + (phi2/phi1) * (dphi1_dxi / dphi2_dxi)) + phi2 * (1 - phi2) * (Pe1 / Pe2)**3)) / ((1-phi1)*((phi1**2) * (phi2 / phi1) * (dphi1_dxi / dphi2_dxi)) + 2 * phi1 * phi2 * ((Pe1 / Pe2)**3) + (phi2**2) * ((Pe1 / Pe2)**6) * (phi1 / phi2) * (dphi2_dxi / dphi1_dxi))
        term2_phi1 = (phi1 + ((Pe1 / Pe2)**3) * phi2) * Z(phi1, phi2)

        term1_phi2 = (K(phi1, phi2) * phi2 * ((phi2 * (1-phi2)) * (phi1/phi2) * (dphi2_dxi / dphi1_dxi) - phi2 * phi1 * (((Pe2 / Pe1)**3) + (phi1/phi2) * (dphi2_dxi / dphi1_dxi)) + phi1 * (1 - phi1) * (Pe2 / Pe1)**3)) / ((1-phi2)*((phi2**2) * (phi1 / phi2) * (dphi2_dxi / dphi1_dxi)) + 2 * phi2 * phi1 * ((Pe2 / Pe1)**3) + (phi1**2) * ((Pe2 / Pe1)**6) * (phi2 / phi1) * (dphi1_dxi / dphi2_dxi))
        term2_phi2 = (phi2 + ((Pe2 / Pe1)**3) * phi1) * Z(phi1, phi2)

        dterm1_phi1_dxi = np.gradient(term1_phi1,xi)
        dterm2_phi1_dxi = np.gradient(term2_phi1,xi)
        
        dterm1_phi2_dxi = np.gradient(term1_phi2,xi)
        dterm2_phi2_dxi = np.gradient(term2_phi2,xi)


        dphi1_dtau = ((1 / (Pe1 * ((1 - t)**2))) * dterm1_phi1_dxi * dterm2_phi1_dxi) - (((xi / (1 - t)) * dphi1_dxi))
        dphi2_dtau = ((1 / (Pe2 * ((1 - t)**2))) * dterm1_phi2_dxi * dterm2_phi2_dxi) - (((xi / (1 - t)) * dphi2_dxi))

        dphi1_dtau[0] = 0
        dphi2_dtau[0] = 0

        return np.concatenate((dphi1_dtau, dphi2_dtau))




    def terminal_condition(t,y,Pe1,Pe2,xi):
        phi1 = y[:N]
        phi2 = y[N:]
        if np.max(phi1) + np.max(phi2) >= 0.64:
            return 0
        return 1


    terminal_condition.terminal = True
    terminal_condition.direction = -1







    y0 = np.concatenate((phi1_0, phi2_0))

    tau_whatch = [0.17,0.34,0.51,0.683]


    sol = solve_ivp(
        pde_system,
        [tau_min, tau_max],
        y0,
        method='RK45',
        t_eval=tau_whatch,
        args=(Pe1,Pe2,xi),
        events=terminal_condition,
        dense_output=True,
        # max_step=1e-3,
        # rtol=1e-6,
        # atol=1e-12
    )

    phi1_sol = sol.y[:N]
    phi2_sol = sol.y[N:]
    print(sol.t)
    print(sol.status)
    print(sol.message)

    num = len(sol.t)


    if num >= 1:
        x  = []
        y1 = []
        y2 = []
        
        for i in range(num):
        # 036  [0.17  0.34  0.51  0.683]
            x.append(xi * (1-tau_whatch[i]))
            y1.append(phi1_sol[:, 0] / (1 - E * tau_whatch[i]))
            y2.append(phi2_sol[:, 0] / (1 - E * tau_whatch[i]))

        plt.figure(figsize=(12, 6))
        for i in range(num):
            # plt.plot(xi, phi1_sol[:, i], label=f'tau={sol.t[i]:.2f}--phi1')
            # plt.plot(xi, phi2_sol[:, i], label=f'tau={sol.t[i]:.2f}--phi2',linestyle='--')
            
            plt.plot(x[i], y1[i], label=f'tau={tau_whatch[i]:.2f}--phi1')
            plt.plot(x[i], y2[i], label=f'tau={tau_whatch[i]:.2f}--phi2',linestyle='--')
            # plt.show()
            
        plt.xlabel(' z/H = xi * (1-tau)')
        plt.ylabel('phi')
        plt.title(f'phi1 & phi2 over z/H  in (Pe1={Pe1}, Pe2={Pe2}, E={E})')
        plt.legend()
        plt.show()
    else:
        print("Solver did not progress beyond the initial state.")









N = 6
E = 1
# phi1_0 = np.linspace(0.101, 0.209, N)
# phi2_0 = np.linspace(0.101, 0.209, N)

phi1_0 = np.array([0.114,0.115,0.118,0.12,0.128,0.138])
phi2_0 = np.array([0.108,0.109,0.111,0.12,0.132,0.152])




# fig 4
grap_phi_1(Pe1=0.18,Pe2=0.37,E=E,phi1_0=phi1_0,phi2_0=phi2_0,N=N)

# fig 5
grap_phi_1(Pe1=0.7,Pe2=0.14,E=E,phi1_0=phi1_0,phi2_0=phi2_0,N=N)

# fig 6
grap_phi_1(Pe1=2.8,Pe2=5.6,E=E,phi1_0=phi1_0,phi2_0=phi2_0,N=N)


Pe1Pe2 = [np.sqrt(0.18 * 0.37) , np.sqrt(0.7 * 1.4),np.sqrt(2.8 * 5.6)]
Pe2_div_Pe1 = [0.37 / 0.18, 1.4 / 0.7, 5.6 / 2.8]




# fig 8 

    # phi2_sol = sol.y[N:]



# fig 9
grap_phi_1(Pe1=0.7,Pe2=0.14,E=E,phi1_0=phi1_0,phi2_0=phi2_0,N=N)

# fig 10
grap_phi_1(Pe1=0.7,Pe2=0.14,E=E,phi1_0=phi1_0,phi2_0=phi2_0,N=N)
    
    
    
    