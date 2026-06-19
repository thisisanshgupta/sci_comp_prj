from sympy import *
from sympy.plotting import plot, plot_parametric

# Projwctile Motion

t = symbols('t', positive=True)
g, k, m = symbols('g k m', positive=True)
v0, theta = symbols('v0 theta', positive=True)

x = Function('x')
z = Function('z')

eq_x = Eq(
    m * diff(x(t), t, 2),
    -k * diff(x(t), t)
)

eq_z = Eq(
    m * diff(z(t), t, 2),
    -m*g - k*diff(z(t), t)
)

print("\nHorizontal Equation:")
pprint(eq_x);

print("\nVertical Equation:")
pprint(eq_z);

ics_x = {
    x(0): 0,
    diff(x(t), t).subs(t, 0): v0*cos(theta)
}

ics_z = {
    z(0): 0,
    diff(z(t), t).subs(t, 0): v0*sin(theta)
}

sol_x = dsolve(eq_x, ics=ics_x)
sol_z = dsolve(eq_z, ics=ics_z)

x_t = simplify(sol_x.rhs)
z_t = simplify(sol_z.rhs)

print("\nSolution for x(t):")
pprint(sol_x);

print("\nSolution for z(t):")
pprint(sol_z);

params = {
    g: 9.81,
    k: 0.15,
    m: 1.0,
    v0: 40,
    theta: pi/4
}

x_num = x_t.subs(params)
z_num = z_t.subs(params)

print("\nSimplified x(t):")
pprint(x_num);

print("\nSimplified z(t):")
pprint(z_num);

eq_x_parabolic = Eq(
    m * diff(x(t), t, 2),
    -0 * diff(x(t), t)
)

eq_z_parabolic = Eq(
    m * diff(z(t), t, 2),
    -m*g - 0 * diff(z(t), t)
)

sol_x_parabolic = dsolve(eq_x_parabolic, ics=ics_x)
sol_z_parabolic = dsolve(eq_z_parabolic, ics=ics_z)

x_t_parabolic = simplify(sol_x_parabolic.rhs)
z_t_parabolic = simplify(sol_z_parabolic.rhs)

params_parabolic_numerical = {
    g: 9.81,
    m: 1.0,
    v0: 40,
    theta: pi/4
}
x_num_parabolic = x_t_parabolic.subs(params_parabolic_numerical)
z_num_parabolic = z_t_parabolic.subs(params_parabolic_numerical)

print("\nSimplified x(t) for parabolic trajectory (k=0):")
pprint(x_num_parabolic);

print("\nSimplified z(t) for parabolic trajectory (k=0):")
pprint(z_num_parabolic);

p1 = plot(
    x_num,
    (t, 0, 10),
    show=False,
    line_color='blue',
    title='Horizontal Position x(t)',
    xlabel='Time',
    ylabel='x(t)'
)

p2 = plot(
    z_num,
    (t, 0, 10),
    show=False,
    line_color='red',
    title='Vertical Position z(t)',
    xlabel='Time',
    ylabel='z(t)'
)

p3 = plot_parametric(
     x_num,
     z_num,
     (t, 0, 10),
     show=False,
     line_color='green',
     title='Trajectory of Particle (with air resistance)',
     xlabel = 'x',
     ylabel = 'z'
)

p4 = plot_parametric(
     x_num_parabolic,
     z_num_parabolic,
     (t, 0, 10),
     show=False,
     line_color='purple',
     title='Parabolic Trajectory (k=0)',
     xlabel = 'x',
     ylabel = 'z'
)

p1.show()
p2.show()
p3.show()
p4.show()


t = symbols('t', positive=True)

g, k, m = symbols('g k m', positive=True)

n = symbols('n', positive=True)

v = Function('v')

eq = Eq(
    m * diff(v(t), t),
    -m*g - k*(v(t))**n
)

eq_linear = eq.subs(n, 1)

ics = {
    v(0): 0
}

sol = dsolve(eq_linear, ics=ics)

v_t = simplify(sol.rhs)

params = {
    g: 9.81,
    k: 1,
    m: 1
}

v_num = v_t.subs(params)

plot(
    v_num,
    (t, 0, 5),
    title="Velocity vs Time with Linear Drag",
    xlabel="Time",
    ylabel="Velocity"
)


t, g, m, c_sym, v0, theta = symbols('t g m c_sym v0 theta', positive=True)
t_dummy = symbols('t_dummy')

v0x = v0 * cos(theta)
v0z = v0 * sin(theta)

gamma_init = 1 / sqrt(1 - (v0**2) / c_sym**2)

px_init = m * v0x * gamma_init
pz_init = m * v0z * gamma_init

px_t = px_init
pz_t = pz_init - m * g * t

Gamma_t = sqrt(1 + (px_t**2 + pz_t**2) / (m**2 * c_sym**2))

vx_t = px_t / (m * Gamma_t)
vz_t = pz_t / (m * Gamma_t)

x_t_rel = integrate(vx_t.subs(t, t_dummy), (t_dummy, 0, t))
z_t_rel = integrate(vz_t.subs(t, t_dummy), (t_dummy, 0, t))

params_rel_toy_c = {
    g: 9.81,
    m: 1.0,
    c_sym: 100,
    v0: 90,
    theta: pi/4
}

x_num_rel = x_t_rel.subs(params_rel_toy_c)
z_num_rel = z_t_rel.subs(params_rel_toy_c)

v0x_classical = params_rel_toy_c[v0] * cos(params_rel_toy_c[theta])
v0z_classical = params_rel_toy_c[v0] * sin(params_rel_toy_c[theta])
g_classical = params_rel_toy_c[g]

x_num_parabolic_compare = v0x_classical * t
z_num_parabolic_compare = v0z_classical * t - (g_classical * t**2) / 2

plot_t_range = (t, 0, 20)

p_rel = plot_parametric(
    x_num_rel,
    z_num_rel,
    plot_t_range,
    show=False,
    line_color='blue',
    label='Relativistic'
)

p_classical = plot_parametric(
    x_num_parabolic_compare,
    z_num_parabolic_compare,
    plot_t_range,
    show=False,
    line_color='red',
    label='Classical'
)

p_rel.extend(p_classical)
p_rel.title = 'Relativistic vs Classical Projectile Motion'
p_rel.xlabel = 'x (m)'
p_rel.ylabel = 'z (m)'
p_rel.legend = True
p_rel.show()
