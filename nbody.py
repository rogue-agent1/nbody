#!/usr/bin/env python3
"""nbody - N-body gravitational simulation with Verlet integration."""
import sys, json, math

class Body:
    def __init__(self, mass, x, y, vx=0, vy=0):
        self.mass = mass; self.x = x; self.y = y; self.vx = vx; self.vy = vy; self.ax = 0; self.ay = 0

def compute_forces(bodies, G=1.0, softening=0.1):
    for b in bodies: b.ax = b.ay = 0
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            dx = bodies[j].x - bodies[i].x; dy = bodies[j].y - bodies[i].y
            r2 = dx*dx + dy*dy + softening*softening
            r = math.sqrt(r2); f = G * bodies[i].mass * bodies[j].mass / r2
            fx, fy = f*dx/r, f*dy/r
            bodies[i].ax += fx/bodies[i].mass; bodies[i].ay += fy/bodies[i].mass
            bodies[j].ax -= fx/bodies[j].mass; bodies[j].ay -= fy/bodies[j].mass

def step(bodies, dt, G=1.0):
    for b in bodies:
        b.x += b.vx*dt + 0.5*b.ax*dt*dt
        b.y += b.vy*dt + 0.5*b.ay*dt*dt
    old_ax = [(b.ax, b.ay) for b in bodies]
    compute_forces(bodies, G)
    for i, b in enumerate(bodies):
        b.vx += 0.5*(old_ax[i][0]+b.ax)*dt
        b.vy += 0.5*(old_ax[i][1]+b.ay)*dt

def total_energy(bodies, G=1.0):
    ke = sum(0.5*b.mass*(b.vx**2+b.vy**2) for b in bodies)
    pe = 0
    for i in range(len(bodies)):
        for j in range(i+1, len(bodies)):
            dx = bodies[j].x-bodies[i].x; dy = bodies[j].y-bodies[i].y
            pe -= G*bodies[i].mass*bodies[j].mass/math.sqrt(dx*dx+dy*dy+0.01)
    return ke + pe

def main():
    bodies = [Body(100,0,0), Body(1,10,0,0,3.16), Body(1,-10,0,0,-3.16), Body(0.5,0,7,2.5,0)]
    compute_forces(bodies)
    e0 = total_energy(bodies)
    print("N-body simulation demo\n")
    for t in range(1000): step(bodies, 0.01)
    e1 = total_energy(bodies)
    print(f"  4 bodies, 1000 steps (dt=0.01)")
    print(f"  Energy: {e0:.4f} -> {e1:.4f} (drift: {abs(e1-e0)/abs(e0)*100:.2f}%)")
    for i, b in enumerate(bodies):
        print(f"  Body {i}: pos=({b.x:.2f},{b.y:.2f}) vel=({b.vx:.2f},{b.vy:.2f})")

if __name__ == "__main__":
    main()
