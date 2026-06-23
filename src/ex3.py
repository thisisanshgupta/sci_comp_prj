# Double Pendulum and Coupled Pendulum (with a spring) using Lagrangian mechanics 

import sympy as sp

# Double Pendulum 
t = sp.symbols('t')

theta1 = sp.Function('theta1')(t)
theta2 = sp.Function('theta2')(t)

m1,m2,l1,l2,g = sp.symbols(
    'm1 m2 l1 l2 g'
)

# coordinates

x1 = l1*sp.sin(theta1)
y1 = -l1*sp.cos(theta1)

x2 = x1 + l2*sp.sin(theta2)
y2 = y1 - l2*sp.cos(theta2)

vx1 = sp.diff(x1,t)
vy1 = sp.diff(y1,t)

vx2 = sp.diff(x2,t)
vy2 = sp.diff(y2,t)

T = (
    sp.Rational(1,2)*m1*(vx1**2 + vy1**2)
    +
    sp.Rational(1,2)*m2*(vx2**2 + vy2**2)
)

V = m1*g*y1 + m2*g*y2

L = T - V

eq1 = (
    sp.diff(
        sp.diff(L,sp.diff(theta1,t)),
        t
    )
    -
    sp.diff(L,theta1)
)

eq2 = (
    sp.diff(
        sp.diff(L,sp.diff(theta2,t)),
        t
    )
    -
    sp.diff(L,theta2)
)

eq1 = sp.simplify(eq1)
eq2 = sp.simplify(eq2)

theta1dd = sp.diff(theta1,t,2)
theta2dd = sp.diff(theta2,t,2)

sol = sp.solve(
    [eq1,eq2],
    [theta1dd,theta2dd]
)

a1 = sol[theta1dd]
a2 = sol[theta2dd]

sp.pprint(a1)
sp.pprint(a2)

# Coupled Pendulum 
t = sp.symbols('t')

theta1 = sp.Function('theta1')(t)
theta2 = sp.Function('theta2')(t)

m,l,g,k = sp.symbols(
    'm l g k'
)

T = (
    sp.Rational(1,2)*m*l**2*sp.diff(theta1,t)**2
    +
    sp.Rational(1,2)*m*l**2*sp.diff(theta2,t)**2
)

V = (
    sp.Rational(1,2)*m*g*l*theta1**2
    +
    sp.Rational(1,2)*m*g*l*theta2**2
    +
    sp.Rational(1,2)*k*l**2
    *(theta2-theta1)**2
)

L = T - V

EL1 = (
    sp.diff(
        sp.diff(L,sp.diff(theta1,t)),
        t
    )
    -
    sp.diff(L,theta1)
)

EL2 = (
    sp.diff(
        sp.diff(L,sp.diff(theta2,t)),
        t
    )
    -
    sp.diff(L,theta2)
)

theta1dd = sp.diff(theta1,t,2)
theta2dd = sp.diff(theta2,t,2)

sol = sp.solve(
    [EL1,EL2],
    [theta1dd,theta2dd]
)

sp.pprint(sol)