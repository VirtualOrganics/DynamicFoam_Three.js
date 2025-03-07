# DynamicFoam Simulation

This is a Three.js implementation of the [DynamicFoam](https://github.com/weigert/DynamicFoam) project by Nicholas McDonald, based on a concept by Michel van de Gaer.

## Description

DynamicFoam creates a "zero-player game" simulation on a foam-like mesh structure derived from a Delaunay triangulation. The system simulates particles flowing along the edges of the foam, which influence the structure of the foam over time, creating interesting patterns similar to reaction-diffusion systems.

## Features

- Interactive 3D visualization using Three.js
- Real-time simulation of foam dynamics
- Three types of triangle center foams: circumcenter, barycenter, and incenter
- Adjustable simulation parameters
- Fully browser-based implementation

## How It Works

1. A set of points are distributed in 2D space
2. Delaunay triangulation is performed on these points
3. Triangle centers are calculated (circumcenter, barycenter, or incenter)
4. Adjacent triangle centers are connected to form the foam structure
5. Particles flow along the edges of the foam, bouncing at the vertices
6. The accumulated flow along edges influences the structure, causing contraction or expansion
7. Over time, patterns emerge in the foam structure

## Usage

### Running the Simulation

Simply open `index.html` in a modern web browser. The simulation will start automatically.

### Controls

- **Mouse**: Drag to rotate, scroll to zoom
- **Space**: Pause/Resume simulation
- **R**: Reset simulation with current parameters

### User Interface

The UI panel on the right allows you to adjust various parameters:

**Simulation:**
- Number of Points: Controls the density of the triangulation
- Center Type: Choose between circumcenter, barycenter, or incenter

**Dynamics:**
- Flow Speed: Controls the speed of particles along edges
- Flow Strength: Determines how much flow accumulates at bounces
- Expansion Threshold: Flow level at which edges expand instead of contract
- Contraction Factor: How much edges contract under normal flow
- Expansion Factor: How much edges expand under high flow
- Equilibrium Distance: Target distance between foam vertices

**Visualization:**
- Show Triangulation: Toggle visibility of the Delaunay triangulation
- Show Foam Edges: Toggle visibility of the foam structure
- Show Flow Particles: Toggle visibility of the flowing particles
- Colors: Customize the colors used in the visualization

## Requirements

- A modern web browser with WebGL support
- No server-side components needed

## Credits

- Original C++ implementation: [weigert/DynamicFoam](https://github.com/weigert/DynamicFoam)
- Concept by Michel van de Gaer
- Delaunay triangulation powered by [Delaunator](https://github.com/mapbox/delaunator)
- 3D visualization using [Three.js](https://threejs.org/)

## License

This project is open source and available under the MIT License. 