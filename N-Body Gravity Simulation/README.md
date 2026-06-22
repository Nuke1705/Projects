# Barnes–Hut N-Body Gravity Simulation

A real-time 2D gravitational N-body simulation written in C++ using OpenGL and ImGui.

This project implements a Barnes–Hut quadtree approximation to accelerate force computation from:

\[
O(N^2) \rightarrow O(N \log N)
\]

allowing simulations with up to 100k particles in real time on CPU.

---

# Features

- Real-time gravitational N-body simulation
- Barnes–Hut quadtree force approximation
- Custom quadtree memory pool allocator
- OpenGL rendering pipeline
- ImGui control/debug interface
- Runtime simulation controls
- Emergent structure formation
- Stable large-scale particle evolution
- Supports 100k+ particles

---

# Simulation Overview

Particles interact through Newtonian gravity:

\[
F \propto \frac{1}{r^2}
\]

To avoid the quadratic complexity of direct pairwise interaction, the simulation uses a Barnes–Hut quadtree structure to approximate distant particle clusters using their center of mass.

This reduces computational complexity significantly while preserving large-scale dynamics.

---

# Emergent Behavior

The simulation exhibits several interesting collective behaviors:

- Filamentary collapse
- Hierarchical clustering
- Violent relaxation
- Approximate virialized equilibrium
- Fractal-like transient structures

Uniform random initial conditions naturally evolve into dense gravitational clusters over time.

---

# Screenshots

## Early Collapse

<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/01a6030a-6de9-403b-828e-ef516717a283" />

## Cluster Formation

<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/715cada6-5cef-4236-9ee2-ff240ee55e5b" />
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/f782a54a-aa5b-42c6-826d-2d1d1bd63653" />
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/7ff8e174-c4a8-4f25-974f-b21ec51684b2" />
<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/28a99c2a-c3f2-4228-b774-07bd3976652a" />


## Stable Virialized State

<img width="1920" height="1080" alt="image" src="https://github.com/user-attachments/assets/49b0c345-22a6-453c-ba83-23f22e2003c1" />

(Drift is due to errors accumulating due to Euler integration which does not conserve energy)
---

# Technologies Used

- C++
- OpenGL
- GLFW
- GLM
- ImGui

---

# Performance

| Particle Count | FPS |
|---|---|
| 10k | ~120 |
| 20k | ~60 |
| 50k | ~30 |
| 100k | ~14 |

Performance scales approximately as:

\[
O(N \log N)
\]

consistent with Barnes–Hut complexity.

---

# Quadtree Structure

Each quadtree node stores:

- Bounding box
- Total mass
- Center of mass
- Child pointers

The tree is rebuilt every frame based on updated particle positions.

---

# Numerical Stability

The simulation includes:

- Gravitational softening
- Minimum node subdivision size
- Recursive force traversal
- Dynamic timestep evolution

Current integration uses Euler integration.

Future improvements include:
- Leapfrog integration
- Symplectic integrators
- Adaptive timesteps
- Proper Barnes–Hut opening criterion

---

# Controls

| Key / UI | Action |
|---|---|
| Pause | Pause simulation |
| Reset | Reset particle distribution |
| ImGui Panel | Runtime diagnostics |

Displayed runtime data:
- Particle count
- FPS
- Quadtree node count
- Timestep

---

# Future Improvements

- Proper Barnes–Hut opening criterion
- GPU compute acceleration
- SIMD optimization
- Parallel tree traversal
- Collision handling
- Galaxy initial conditions
- Energy diagnostics
- Center-of-mass tracking
- Velocity Verlet / Leapfrog integration

---

# Building

## Requirements

- C++17
- OpenGL
- GLFW
- GLM (Please include the glm files inside a folder named glm under the resources folder)
- ImGui
