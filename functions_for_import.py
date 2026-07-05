def eulerStep(h, derivative, x, y):
    return (x + h, y + h*derivative(x,y))

def eulersMethod(states, steps, h, derivative, x0, y0,):
    for _ in range(0, steps):
        x0, y0 = eulerStep(h, derivative, x0, y0)
        states.append((x0,y0))
    return states

def leapfrogStep(h, derivative, x0, y0, x1, y1):
    # Apply leapfrog method
    x2 = x1 + h
    y2 = y0 + 2*h*derivative(x1,y1)
    return x1, y1, x2, y2

def leapfrogMethod(states, steps, h, derivative, x0, y0):
    # Firstly apply one step of Euler's
    x1, y1 = eulerStep(h, derivative, x0, y0)
    states.append((x1,y1))
    # Then apply the leapfrog method for the remaining steps
    for _ in range(0,steps - 1):
        x0, y0, x1, y1 = leapfrogStep(h, derivative, x0, y0, x1, y1)
        states.append((x1,y1))
    return states

def midpointStep(h, derivative, x0, y0):
    # Apply midpoint / RK2 method
    x1 = x0 + h
    xmid = x0 + h/2
    ymid = y0 + h/2 * derivative(x0, y0)
    y1 = y0 + h*derivative(xmid, ymid)
    return (x1, y1)

def midpointMethod(states, steps, h, derivative, x0, y0):
    for _ in range(0, steps):
        x0, y0 = midpointStep(h, derivative, x0, y0)
        states.append((x0,y0))
    return states

def rungeKuttaStep(h, derivative, x0, y0):
    # Apply Runge-Kutta / RK4 method
    k1 = derivative(x0, y0)
    k2 = derivative(x0 + h/2, y0 + h/2 * k1)
    k3 = derivative(x0 + h/2, y0 + h/2 * k2)
    k4 = derivative(x0 + h, y0 + h * k3)
    x1 = x0 + h
    y1 = y0 + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    return (x1, y1)

def rungeKuttaMethod(states, steps, h, derivative, x0, y0):
    for _ in range(0, steps):
        x0, y0 = rungeKuttaStep(h, derivative, x0, y0)
        states.append((x0,y0))
    return states

def runAllMethods(initialiseStates, states, steps, h, derivative, x0, y0):
    eulerStates = eulersMethod(states, steps, h, derivative, x0, y0)
    x0, y0, h, states, steps = initialiseStates()
    leapfrogStates = leapfrogMethod(states, steps, h, derivative, x0, y0)
    x0, y0, h, states, steps = initialiseStates()
    midpointStates = midpointMethod(states, steps, h, derivative, x0, y0)
    x0, y0, h, states, steps = initialiseStates()
    rungeKuttaStates = rungeKuttaMethod(states, steps, h, derivative, x0, y0)
    return eulerStates, leapfrogStates, midpointStates, rungeKuttaStates

def compareMethods(eulerStates, leapfrogStates, midpointStates, rungeKuttaStates):
    print("     x     |  y (Euler)   | y (Leapfrog) | y (Midpoint) |  y (RK4)   ")
    print("-----------+--------------+--------------+--------------+-------------")
    for i in range(len(eulerStates)):
        x = eulerStates[i][0]
        y_euler = eulerStates[i][1]
        y_leapfrog = leapfrogStates[i][1]
        y_midpoint = midpointStates[i][1]
        y_rk4 = rungeKuttaStates[i][1]
        print(f"{x:10.4f} | {y_euler:12.6f} | {y_leapfrog:12.6f} | {y_midpoint:12.6f} | {y_rk4:12.6f}")