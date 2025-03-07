import { GUI } from 'three/addons/libs/lil-gui.module.min.js';

/**
 * UIController class handles the user interface for adjusting simulation parameters
 */
class UIController {
    constructor(container, foam, renderer) {
        this.container = container;
        this.foam = foam;
        this.renderer = renderer;
        
        // GUI parameters
        this.params = {
            // Simulation parameters
            numPoints: foam.numPoints,
            centerType: foam.centerType,
            flowSpeed: foam.flowSpeed,
            flowStrength: foam.flowStrength,
            expansionThreshold: foam.expansionThreshold,
            contractionFactor: foam.contractionFactor,
            expansionFactor: foam.expansionFactor,
            equilibriumDistance: foam.equilibriumDistance,
            
            // Visualization parameters
            showDelaunayEdges: true,     // The backbone of our model
            showFoamEdges: true,         // Voronoi edges
            showFlowParticles: true,     // Flow particles
            
            // Actions
            isRunning: true,
            toggleSimulation: () => {
                this.params.isRunning = window.app.toggleSimulation();
            },
            resetSimulation: () => this.resetSimulation()
        };
        
        // Initialize GUI
        this.gui = new GUI({ container: this.container });
        this.setupGUI();
    }
    
    /**
     * Set up the GUI controls
     */
    setupGUI() {
        // Create folders
        const simulationFolder = this.gui.addFolder('Simulation');
        const centerTypeFolder = simulationFolder.addFolder('Center Type');
        const dynamicsFolder = simulationFolder.addFolder('Dynamics');
        const particleSettingsFolder = this.gui.addFolder('Particle Settings');
        const visualizationFolder = this.gui.addFolder('Visualization');
        const actionsFolder = this.gui.addFolder('Actions');
        
        // Center type selection
        centerTypeFolder.add(this.params, 'centerType', ['circumcenter', 'barycenter', 'incenter'])
            .name('Center Type')
            .onChange(value => {
                this.foam.setCenterType(value);
                this.updateRenderer();
            });
        
        // Simulation parameters
        simulationFolder.add(this.params, 'numPoints', 50, 500, 10)
            .name('Number of Points');
        
        // Dynamics parameters
        dynamicsFolder.add(this.params, 'flowSpeed', 0.1, 2.0, 0.1)
            .name('Flow Speed')
            .onChange(value => {
                this.foam.setDynamicsParams({ flowSpeed: value });
            });
        
        dynamicsFolder.add(this.params, 'flowStrength', 0.001, 0.1, 0.001)
            .name('Flow Strength')
            .onChange(value => {
                this.foam.setDynamicsParams({ flowStrength: value });
            });
        
        dynamicsFolder.add(this.params, 'expansionThreshold', 0.1, 1.0, 0.1)
            .name('Expansion Threshold')
            .onChange(value => {
                this.foam.setDynamicsParams({ expansionThreshold: value });
            });
        
        dynamicsFolder.add(this.params, 'contractionFactor', 0.8, 0.99, 0.01)
            .name('Contraction Factor')
            .onChange(value => {
                this.foam.setDynamicsParams({ contractionFactor: value });
            });
        
        dynamicsFolder.add(this.params, 'expansionFactor', 1.01, 1.2, 0.01)
            .name('Expansion Factor')
            .onChange(value => {
                this.foam.setDynamicsParams({ expansionFactor: value });
            });
        
        dynamicsFolder.add(this.params, 'equilibriumDistance', 10, 100, 5)
            .name('Equilibrium Distance')
            .onChange(value => {
                this.foam.setDynamicsParams({ equilibriumDistance: value });
            });
            
        // Particle settings
        particleSettingsFolder.add(this.foam, 'particleReleaseInterval', 0.1, 2.0, 0.1)
            .name('Release Interval (s)');
            
        particleSettingsFolder.add(this.foam, 'particleLifetime', 0.5, 5.0, 0.5)
            .name('Particle Lifetime (s)');
            
        particleSettingsFolder.add(this.foam, 'maxAngleForTraversal', 0, Math.PI, 0.1)
            .name('Max Angle (rad)')
            .onChange(value => {
                // Clamp between 0 and Pi (0-180 degrees)
                this.foam.maxAngleForTraversal = Math.max(0, Math.min(Math.PI, value));
            });
        
        // Visualization parameters
        visualizationFolder.add(this.params, 'showDelaunayEdges')
            .name('Show Delaunay Edges')
            .onChange(value => {
                // Only show Delaunay edges, not triangulation
                this.renderer.setVisibility({ showDelaunayEdges: value });
            });
        
        visualizationFolder.add(this.params, 'showFoamEdges')
            .name('Show Voronoi Mesh')
            .onChange(value => {
                this.renderer.setVisibility({ showFoamEdges: value });
            });
        
        visualizationFolder.add(this.params, 'showFlowParticles')
            .name('Show Flow Particles')
            .onChange(value => {
                this.renderer.setVisibility({ showFlowParticles: value });
            });
        
        // Color controls
        const colorParams = {
            backgroundColor: '#111111',
            delaunayEdgesColor: '#7a7a7a',
            voronoiMeshColor: '#2288ff',
            flowParticlesColor: '#ff8822'
        };
        
        const colorFolder = visualizationFolder.addFolder('Colors');
        
        colorFolder.addColor(colorParams, 'backgroundColor')
            .name('Background')
            .onChange(value => {
                this.renderer.setColors({ background: value });
            });
        
        colorFolder.addColor(colorParams, 'delaunayEdgesColor')
            .name('Delaunay Edges')
            .onChange(value => {
                this.renderer.setColors({ delaunayLines: value });
            });
        
        colorFolder.addColor(colorParams, 'voronoiMeshColor')
            .name('Voronoi Mesh')
            .onChange(value => {
                this.renderer.setColors({ foamLines: value });
            });
        
        colorFolder.addColor(colorParams, 'flowParticlesColor')
            .name('Flow Particles')
            .onChange(value => {
                this.renderer.setColors({ flowParticles: value });
            });
        
        // Actions
        actionsFolder.add(this.params, 'toggleSimulation')
            .name('Pause/Play Simulation');
        actionsFolder.add(this.params, 'resetSimulation').name('Reset Simulation');
        
        // Open folders by default
        simulationFolder.open();
        dynamicsFolder.open();
        particleSettingsFolder.open();
        visualizationFolder.open();
        actionsFolder.open();
    }
    
    /**
     * Reset the simulation with current parameters
     */
    resetSimulation() {
        this.foam.reset(this.params.numPoints);
        this.updateRenderer();
    }
    
    /**
     * Update the renderer with new foam data
     */
    updateRenderer() {
        const foamData = this.foam.getGeometryData();
        this.renderer.updateFoamData(foamData);
    }
    
    /**
     * Update the UI with current foam parameters
     */
    updateUI() {
        // This can be used to sync the UI if foam parameters change from outside
        this.params.numPoints = this.foam.numPoints;
        this.params.centerType = this.foam.centerType;
        this.params.flowSpeed = this.foam.flowSpeed;
        this.params.flowStrength = this.foam.flowStrength;
        this.params.expansionThreshold = this.foam.expansionThreshold;
        this.params.contractionFactor = this.foam.contractionFactor;
        this.params.expansionFactor = this.foam.expansionFactor;
        this.params.equilibriumDistance = this.foam.equilibriumDistance;
        
        // Update GUI display
        this.gui.controllers.forEach(controller => controller.updateDisplay());
    }
}

export default UIController; 