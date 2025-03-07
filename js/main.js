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
        
        // Initial render
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
     * Update flow particle visualization
     */
    updateFlowParticles() {
        const foamData = this.foam.getGeometryData();
        this.renderer.updateFlowParticles(
            foamData.centers,
            foamData.edges,
            foamData.flows
        );
    }
    
    /**
     * Update foam edges visualization
     */
    updateFoamEdges() {
        const foamData = this.foam.getGeometryData();
        this.renderer.updateFoamEdges(
            foamData.centers,
            foamData.edges
        );
    }
    
    /**
     * Animation loop
     */
    animate(time) {
        requestAnimationFrame(this.animate.bind(this));
        
        if (!this.isRunning) return;
        
        // Calculate time delta
        if (!time) time = 0;
        const deltaTime = (time - this.lastTime) / 1000; // in seconds
        this.lastTime = time;
        
        // Skip if deltaTime is too large (e.g., after tab switch)
        if (deltaTime > 0.1) return;
        
        // Update foam simulation
        this.foam.update(deltaTime);
        
        // Update visualizations
        this.updateFlowParticles();
        this.updateFoamEdges();
        
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
