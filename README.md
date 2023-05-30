# Modélisation d'un double pendule

Projet ayant pour objectif de modéliser un double pendule à l'aide de MatLab.

Double pendule : 

![image](https://user-images.githubusercontent.com/58084848/177021215-03486d6f-de89-4f58-8ff9-d54e603a53e3.png)

Modélisation du mouvement selon plusieurs méthodes numériques :

![image](https://user-images.githubusercontent.com/58084848/177039991-a71f5317-984d-446c-b398-3745c9f84eac.png)

Le mouvement d'un double pendule peut être modélisé à l'aide d'équations différentielles couplées qui décrivent l'évolution des angles et des vitesses angulaires des deux pendules au fil du temps.

Les équations qui décrivent le mouvement d'un double pendule sont non-linéaires et peuvent être difficiles à résoudre analytiquement. Cependant, elles peuvent être résolues numériquement à l'aide de méthodes numériques.

Voici les équations qui décrivent le mouvement d'un double pendule :

```
θ1'' = (-g*(2m1+m2)sin(θ1) - m2gsin(θ1-2θ2) - 2sin(θ1-θ2)m2(θ2'^2L2+θ1'^2L1cos(θ1-θ2))) / (L1(2m1+m2-m2cos(2θ1-2θ2)))

θ2'' = (2sin(θ1-θ2)(θ1'^2L1(m1+m2)+g*(m1+m2)cos(θ1)+θ2'^2L2m2cos(θ1-θ2))) / (L2*(2m1+m2-m2cos(2θ1-2θ2)))
```

où :
- θ1 et θ2 sont les angles des deux pendules par rapport à la verticale.
- m1 et m2 sont les masses des deux pendules.
- L1 et L2 sont les longueurs des deux pendules.
- g est l'accélération due à la gravité.
- θ1' et θ2' sont les vitesses angulaires des deux pendules.
- θ1'' et θ2'' sont les accélérations angulaires des deux pendules.

En résolvant ces équations numériquement, on peut obtenir une représentation du mouvement du double pendule dans le temps. Le mouvement du double pendule est souvent chaotique et difficile à prédire à l'avance, ce qui en fait l'exemple classique de système dynamique non-linéaire.

Plus d'explications dans le rapport "PROJMM8.pdf".  

## Exécution

Lancez le programme "PROJMM8.m" dans MatLab et décommentez les graphes que vous voulez voir.

Le programme python est un exemple de code qui utilise la méthode de Runge-Kutta pour résoudre les équations différentielles du double pendule et tracer la trajectoire des deux pendules. Ce code utilise la bibliothèque scipy pour résoudre les équations différentielles.

``` python
# Fonction qui calcule les dérivées de theta1, theta2, theta1', et theta2' à partir des valeurs actuelles
def derivatives(y, t):
    theta1, theta2, omega1, omega2 = y
    delta_theta = theta2 - theta1
    delta_omega = omega2 - omega1
    
    # Calcul des dérivées
    dtheta1_dt = omega1
    dtheta2_dt = omega2
    domega1_dt = (m2*L1*omega1**2*np.sin(delta_theta)*np.cos(delta_theta) + m2*g*np.sin(theta2)*np.cos(delta_theta) + m2*L2*omega2**2*np.sin(delta_theta) - (m1+m2)*g*np.sin(theta1)) / ((m1+m2)*L1 - m2*L1*np.cos(delta_theta)**2)
    domega2_dt = (-m2*L2*omega2**2*np.sin(delta_theta)*np.cos(delta_theta) + (m1+m2)*(g*np.sin(theta1)*np.cos(delta_theta) - L1*omega1**2*np.sin(delta_theta) - g*np.sin(theta2)) / L2) / ((m1+m2)*L2 - m2*L2*np.cos(delta_theta)**2)
    
    return dtheta1_dt, dtheta2_dt, domega1_dt, domega2_dt
 ```

![image](https://user-images.githubusercontent.com/58084848/234627242-a8f9c598-e79a-4719-a665-b6530bc415da.png)

## Crédits

- [Clément DILLY](https://youtu.be/dQw4w9WgXcQ)
- [Guillaume DUMAS](https://github.com/Nausicaa68)
- [Jérémy GRELAUD](https://github.com/jeremyGrelaud) 
- [Julien HASSOUN](https://youtu.be/dQw4w9WgXcQ)
- [Amaury ROSSIGNOL](https://github.com/Mushurisen)
