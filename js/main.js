import * as THREE from 'three';
import DelaunayFoam from './delaunayFoam.js';
import FoamRenderer from './foamRenderer.js';
import UIController from './uiController.js';

/**
 * Main application class for the DynamicFoam simulation
 */
class DynamicFoamApp {
    constructor() {
        // Set up foam simulation
        this.foam = new DelaunayFoam(800, 600, 150);
        
        // Set up Three.js renderer
        const container = document.getElementById('canvas');
        this.renderer = new FoamRenderer(container, window.innerWidth, window.innerHeight);
        
        // Set up UI controller
        const uiContainer = document.getElementById('gui-container');
        this.ui = new UIController(uiContainer, this.foam, this.renderer);
        
        // Animation properties
        this.lastTime = 0;
        this.isRunning = true;
        this.lastDebugTime = 0;
        this.lastLogTime = 0;
        
        // Initial render
        this.updateRenderer();
        
        // Ensure initial particles are created
        this.foam.initializeFlows();
        this.updateRenderer();
        
        // Start animation loop
        this.animate();
    }
    
    /**
     * Update the renderer with current foam data
     */
    updateRenderer() {
        const foamData = this.foam.getGeometryData();
        this.renderer.updateFoamData(foamData);
    }
    
    /**
     * Update all visualization components
     */
    updateVisualizations() {
        const foamData = this.foam.getGeometryData();
        
        // Update flow particles
        this.renderer.updateFlowParticles(
            foamData.centers,
            foamData.edges,
            foamData.flows
        );
        
        // Update foam edges (Voronoi)
        this.renderer.updateFoamEdges(
            foamData.centers,
            foamData.edges
        );
        
        // Update Delaunay edges
        this.renderer.updateDelaunayEdges(
            foamData.points,
            foamData.delaunayEdges
        );
    }
    
    /**
     * Animation loop
     */
    animate(time) {
        requestAnimationFrame(this.animate.bind(this));
        
        if (!this.isRunning) {
            this.renderer.render();
            return;
        }
        
        // Calculate time delta
        if (!time) time = 0;
        const deltaTime = Math.min((time - this.lastTime) / 1000, 0.05); // Cap at 50ms to avoid large jumps
        this.lastTime = time;
        
        // Skip if deltaTime is too small
        if (deltaTime < 0.001) {
            this.renderer.render();
            return;
        }
        
        // Update foam simulation
        this.foam.update(deltaTime);
        
        // Debug logging every 2 seconds
        if (Math.floor(time / 2000) !== Math.floor(this.lastLogTime / 2000)) {
            const foamData = this.foam.getGeometryData();
            console.log(`Active particles: ${foamData.flows.length}`);
            this.lastLogTime = time;
        }
        
        // Update visualizations
        this.updateVisualizations();
        
        // Render scene
        this.renderer.render();
    }
    
    /**
     * Toggle the simulation
     */
    toggleSimulation() {
        this.isRunning = !this.isRunning;
        return this.isRunning;
    }
    
    /**
     * Reset the simulation
     */
    resetSimulation() {
        this.ui.resetSimulation();
    }
}

// Initialize the application when the page is loaded
window.addEventListener('DOMContentLoaded', () => {
    window.app = new DynamicFoamApp();
    
    // Add keyboard controls
    window.addEventListener('keydown', (event) => {
        if (event.key === ' ') {
            // Space to toggle simulation
            const isRunning = window.app.toggleSimulation();
            console.log(`Simulation ${isRunning ? 'resumed' : 'paused'}`);
        } else if (event.key === 'r') {
            // R to reset
            window.app.resetSimulation();
            console.log('Simulation reset');
        }
    });
    
    console.log("DynamicFoam simulation started");
}); 