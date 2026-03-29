#!/usr/bin/env python3
"""nbody - N-body gravitational simulation."""
import sys, math

class Body:
    def __init__(self, mass, x, y, vx=0, vy=0):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

def simulate(bodies, dt, steps, G=1.0, softening=0.01):
    history = [[(b.x, b.y) for b in bodies]]
    for _ in range(steps):
        # compute forces
        forces = [(0.0, 0.0)] * len(bodies)
        for i in range(len(bodies)):
            fx, fy = 0.0, 0.0
            for j in range(len(bodies)):
                if i == j:
                    continue
                dx = bodies[j].x - bodies[i].x
                dy = bodies[j].y - bodies[i].y
                r = math.sqrt(dx*dx + dy*dy + softening**2)
                f = G * bodies[i].mass * bodies[j].mass / (r * r * r)
                fx += f * dx
                fy += f * dy
            forces[i] = (fx, fy)
        # update velocities and positions
        for i, b in enumerate(bodies):
            b.vx += forces[i][0] / b.mass * dt
            b.vy += forces[i][1] / b.mass * dt
            b.x += b.vx * dt
            b.y += b.vy * dt
        history.append([(b.x, b.y) for b in bodies])
    return history

def total_energy(bodies, G=1.0, softening=0.01):
    ke = sum(0.5 * b.mass * (b.vx**2 + b.vy**2) for b in bodies)
    pe = 0.0
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            dx = bodies[j].x - bodies[i].x
            dy = bodies[j].y - bodies[i].y
            r = math.sqrt(dx*dx + dy*dy + softening**2)
            pe -= G * bodies[i].mass * bodies[j].mass / r
    return ke + pe

def test():
    # two bodies orbiting
    b1 = Body(1000, 0, 0)
    b2 = Body(1, 10, 0, 0, 10)
    bodies = [b1, b2]
    e0 = total_energy(bodies)
    history = simulate(bodies, 0.001, 1000)
    e1 = total_energy(bodies)
    # energy should be roughly conserved (small dt)
    assert abs(e1 - e0) / abs(e0) < 0.05
    assert len(history) == 1001
    # body moved
    assert history[-1][1] != history[0][1]
    print("OK: nbody")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        test()
    else:
        print("Usage: nbody.py test")
