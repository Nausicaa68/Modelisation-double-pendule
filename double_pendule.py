from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Paramètres du système
m1 = 1.0  # masse du premier pendule
m2 = 1.0  # masse du deuxième pendule
L1 = 1.0  # longueur du premier pendule
L2 = 1.0  # longueur du deuxième pendule
g = 9.81  # accélération due à la gravité

# Fonction qui calcule les dérivées de theta1, theta2, theta1', et theta2' à partir des valeurs actuelles


def derivatives(y, t):
    print("Processing ...")
    theta1, theta2, omega1, omega2 = y
    delta_theta = theta2 - theta1
    delta_omega = omega2 - omega1

    # Calcul des dérivées
    dtheta1_dt = omega1
    dtheta2_dt = omega2
    domega1_dt = (m2*L1*omega1**2*np.sin(delta_theta)*np.cos(delta_theta) + m2*g*np.sin(theta2)*np.cos(delta_theta) +
                  m2*L2*omega2**2*np.sin(delta_theta) - (m1+m2)*g*np.sin(theta1)) / ((m1+m2)*L1 - m2*L1*np.cos(delta_theta)**2)
    domega2_dt = (-m2*L2*omega2**2*np.sin(delta_theta)*np.cos(delta_theta) + (m1+m2)*(g*np.sin(theta1)*np.cos(delta_theta) -
                  L1*omega1**2*np.sin(delta_theta) - g*np.sin(theta2)) / L2) / ((m1+m2)*L2 - m2*L2*np.cos(delta_theta)**2)

    print("Calcul des dérivées terminés")
    return dtheta1_dt, dtheta2_dt, domega1_dt, domega2_dt


# Intervalle de temps et conditions initiales
t = np.arange(0, 10, 0.01)
# valeurs initiales pour theta1, theta2, theta1', et theta2'
y0 = np.array([np.pi/2, np.pi/2, 0, 0])

# Résolution des équations différentielles à l'aide de la méthode de Runge-Kutta
sol = odeint(derivatives, y0, t)

# Extraction des valeurs de theta1, theta2, x1, y1, x2, et y2 à partir de la solution
theta1 = sol[:, 0]
theta2 = sol[:, 1]
x1 = L1*np.sin(theta1)
y1 = -L1*np.cos(theta1)
x2 = x1 + L2*np.sin(theta2)
y2 = y1 - L2*np.cos(theta2)

# Tracé de la trajectoire des deux pendules
print("Tracé des trajectoires")
fig, ax = plt.subplots()
ax.plot(x1, y1, 'b', label='Pendule 1')
ax.plot(x2, y2, 'g', label='Pendule 2')
ax.set_xlabel('Position horizontale')
ax.set_ylabel('Position verticale')
ax.set_title('Trajectoire du double pendule')
ax.legend()
plt.show()
