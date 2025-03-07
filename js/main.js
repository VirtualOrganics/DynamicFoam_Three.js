import * as THREE from 'three';
import DelaunayFoam from './delaunayFoam.js';
import FoamRenderer from './foamRenderer.js';
import UIController from './uiController.js';

/**
 * DynamicFoamApp - Main application class
 */
class DynamicFoamApp {
    constructor() {
        this.width = window.innerWidth;
        this.height = window.innerHeight;
        this.container = document.getElementById('container');
        
        // Create foam simulation
        this.foam = new DelaunayFoam(this.width, this.height, 100);
        
        // Create renderer
        this.renderer = new FoamRenderer(this.container, this.width, this.height);
        
        // Create UI controller
        this.uiController = new UIController(this.container, this.foam, this.renderer);
        
        // Animation state
        this.isRunning = true;
        this.lastTime = performance.now() / 1000;
        this.lastLogTime = 0;
        
        // Start animation
        this.animate();
        
        // Initialize foam data to ensure first frame has particles
        this.updateRenderer();
        
        // Ensure particles are created at startup
        this.foam.initializeFlows();
        const foamData = this.foam.getGeometryData();
        this.renderer.updateFoamData(foamData);
        
        console.log("DynamicFoamApp initialized!");
    }
    
    /**
     * Update renderer with current foam data
     */
    updateRenderer() {
        const foamData = this.foam.getGeometryData();
        this.renderer.updateFoamData(foamData);
    }
    
    /**
     * Update visualization options
     */
    updateVisualizations() {
        // Get updated geometry data
        const foamData = this.foam.getGeometryData();
        
        // Update renderer with new data
        this.renderer.updateFoamData(foamData);
    }
    
    /**
     * Animation loop
     */
    animate() {
        requestAnimationFrame(this.animate.bind(this));
        
        if (!this.isRunning) {
            this.renderer.render();
            return;
        }
        
        // Calculate delta time (in seconds)
        const currentTime = performance.now() / 1000;
        const deltaTime = Math.min(0.1, currentTime - this.lastTime); // Cap at 100ms to avoid large jumps
        this.lastTime = currentTime;
        
        // Update the simulation
        this.foam.update(deltaTime);
        
        // Get updated geometry data
        const foamData = this.foam.getGeometryData();
        
        // Debug log flow data
        if (Math.floor(currentTime) % 5 === 0 && Math.floor(currentTime) !== this.lastLogTime) {
            this.lastLogTime = Math.floor(currentTime);
            console.log(`Active particles: ${foamData.flows.length}`);
        }
        
        // Update the renderer with new data
        this.renderer.updateFoamData(foamData);
        
        // Render the scene
        this.renderer.render();
    }
    
    /**
     * Toggle simulation running state
     */
    toggleSimulation() {
        this.isRunning = !this.isRunning;
        if (this.isRunning) {
            this.lastTime = performance.now() / 1000;
        }
    }
    
    /**
     * Reset simulation
     */
    resetSimulation() {
        this.foam.reset();
        this.foam.initializeFlows(); // Ensure particles are immediately created
        this.updateRenderer();
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