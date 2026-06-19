from sympy import *
from sympy.plotting import plot, plot_parametric

# Wave Mechanics
t = symbols('t')

a, b = symbols('a b', positive=True)

w1, w2 = symbols('w1 w2', positive=True)

phi = symbols('phi')

x = a*sin(w1*t + phi)
y = b*sin(w2*t)

params = {
    a: 1,
    b: 1,
    w1: 3,
    w2: 2,
    phi: pi/2
}

x_num = x.subs(params)
y_num = y.subs(params)

plot_parametric(
    (x_num, y_num),
    (t, 0, 2*pi),
    title="Lissajous Figure",
    xlabel="x",
    ylabel="y"
)


x, t = symbols('x t')

v1, v2, v3 = symbols('v1 v2 v3')

k1, k2, k3 = symbols('k1 k2 k3')

y1 = 3*sin(2*pi*v1*t - k1*x)

y2 = 3*sin(2*pi*v2*t - k2*x)

y3 = 3*sin(2*pi*v3*t - k3*x)

y = simplify(y1 + y2 + y3)

params = {
    v1: 9.5,
    v2: 10,
    v3: 10.5,
    k1: 15,
    k2: 15.5,
    k3: 16,
    t: 5
}

y1_num = y1.subs(params)
y2_num = y2.subs(params)
y3_num = y3.subs(params)

y_num = y.subs(params)

p1 = plot(
    y1_num,
    (x, 0, 15),
    show=False,
    title="Wave y1"
)

p2 = plot(
    y2_num,
    (x, 0, 10),
    show=False,
    title="Wave y2"
)

p3 = plot(
    y3_num,
    (x, 0, 10),
    show=False,
    title="Wave y3"
)

p4 = plot(
    y_num,
    (x, 0, 10),
    show=False,
    title="Superposed Wave"
)

p1.show()
p2.show()
p3.show()
p4.show()



