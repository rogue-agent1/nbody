import argparse, math, random

class Body:
    def __init__(self, x, y, vx, vy, mass):
        self.x, self.y = x, y
        self.vx, self.vy = vx, vy
        self.mass = mass

def simulate(bodies, dt, steps, G=6.674e-11):
    history = []
    for step in range(steps):
        if step % (steps//5) == 0:
            history.append([(b.x, b.y) for b in bodies])
        for i, a in enumerate(bodies):
            fx, fy = 0, 0
            for j, b in enumerate(bodies):
                if i == j: continue
                dx, dy = b.x - a.x, b.y - a.y
                r = math.sqrt(dx*dx + dy*dy + 1e-10)
                f = G * a.mass * b.mass / (r * r)
                fx += f * dx / r
                fy += f * dy / r
            a.vx += fx / a.mass * dt
            a.vy += fy / a.mass * dt
        for b in bodies:
            b.x += b.vx * dt
            b.y += b.vy * dt
    return history

def main():
    p = argparse.ArgumentParser(description="N-body simulation")
    p.add_argument("-n", "--bodies", type=int, default=5)
    p.add_argument("-s", "--steps", type=int, default=1000)
    p.add_argument("--solar", action="store_true", help="Solar system demo")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    if args.solar:
        AU = 1.496e11
        bodies = [
            Body(0, 0, 0, 0, 1.989e30),  # Sun
            Body(0.387*AU, 0, 0, 47400, 3.285e23),  # Mercury
            Body(0.723*AU, 0, 0, 35000, 4.867e24),  # Venus
            Body(AU, 0, 0, 29800, 5.972e24),  # Earth
            Body(1.524*AU, 0, 0, 24100, 6.39e23),  # Mars
        ]
        dt = 86400  # 1 day
    else:
        random.seed(args.seed)
        bodies = [Body(random.uniform(-1e10, 1e10), random.uniform(-1e10, 1e10),
                       random.uniform(-1e3, 1e3), random.uniform(-1e3, 1e3),
                       random.uniform(1e24, 1e30)) for _ in range(args.bodies)]
        dt = 3600
    history = simulate(bodies, dt, args.steps)
    for i, frame in enumerate(history):
        print(f"Frame {i}:")
        for j, (x, y) in enumerate(frame):
            name = ["Sun","Mercury","Venus","Earth","Mars"][j] if args.solar and j < 5 else f"Body {j}"
            print(f"  {name}: ({x:.3e}, {y:.3e})")

if __name__ == "__main__":
    main()
