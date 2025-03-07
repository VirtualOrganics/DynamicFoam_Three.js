import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

/**
 * FoamRenderer class handles the visualization of the foam structure using Three.js
 */
class FoamRenderer {
    constructor(container, width, height) {
        this.container = container;
        this.width = width;
        this.height = height;
        
        // Three.js objects
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        
        // Visualization objects
        this.triangles = null;
        this.foamLines = null;
        this.flowParticles = null;
        
        // Materials
        this.triangleMaterial = null;
        this.lineMaterial = null;
        this.particleMaterial = null;
        
        // Colors
        this.backgroundColor = new THREE.Color(0x111111);
        this.triangleColor = new THREE.Color(0x222222);
        this.foamLineColor = new THREE.Color(0x2288ff);
        this.flowParticleColor = new THREE.Color(0xff8822);
        
        // Initialize
        this.init();
    }
    
    /**
     * Initialize the Three.js scene, camera, and renderer
     */
    init() {
        // Create scene
        this.scene = new THREE.Scene();
        this.scene.background = this.backgroundColor;
        
        // Create camera
        this.camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 2000);
        this.camera.position.z = 800;
        
        // Create renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);
        
        // Add orbit controls
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        
        // Create materials
        this.triangleMaterial = new THREE.MeshBasicMaterial({
            color: 0x444444,  // Lighter gray for more contrast
            wireframe: true,
            transparent: true,
            opacity: 0.4      // Slightly more opaque
        });
        
        this.lineMaterial = new THREE.LineBasicMaterial({
            color: this.foamLineColor,
            linewidth: 1
        });
        
        this.particleMaterial = new THREE.PointsMaterial({
            color: this.flowParticleColor,
            size: 3,
            sizeAttenuation: true
        });
        
        // Handle window resize
        window.addEventListener('resize', this.onWindowResize.bind(this));
    }
    
    /**
     * Resize handler
     */
    onWindowResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }
    
    /**
     * Update scene with new foam data
     */
    updateFoamData(foamData) {
        // Remove previous objects
        if (this.triangles) this.scene.remove(this.triangles);
        if (this.foamLines) this.scene.remove(this.foamLines);
        if (this.flowParticles) this.scene.remove(this.flowParticles);
        
        // Create triangulation visualization
        this.createTriangulation(foamData.points, foamData.triangles);
        
        // Create foam edges
        this.createFoamEdges(foamData.centers, foamData.edges);
        
        // Create flow particles
        this.createFlowParticles(foamData.centers, foamData.edges, foamData.flows);
    }
    
    /**
     * Create visualization of the Delaunay triangulation
     */
    createTriangulation(points, triangles) {
        const geometry = new THREE.BufferGeometry();
        
        // Convert points to Vector3 array for Three.js
        const vertices = [];
        for (let i = 0; i < points.length; i += 2) {
            vertices.push(points[i], points[i + 1], 0);
        }
        
        geometry.setIndex(triangles);
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        
        this.triangles = new THREE.Mesh(geometry, this.triangleMaterial);
        this.scene.add(this.triangles);
    }
    
    /**
     * Create visualization of foam edges
     */
    createFoamEdges(centers, edges) {
        const geometry = new THREE.BufferGeometry();
        
        // Create line segments from edges
        const vertices = [];
        for (const edge of edges) {
            const from = centers[edge.from];
            const to = centers[edge.to];
            
            vertices.push(from[0], from[1], 0);
            vertices.push(to[0], to[1], 0);
        }
        
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        
        this.foamLines = new THREE.LineSegments(geometry, this.lineMaterial);
        this.scene.add(this.foamLines);
    }
    
    /**
     * Create visualization of flow particles
     */
    createFlowParticles(centers, edges, flows) {
        const geometry = new THREE.BufferGeometry();
        
        // Create particles at flow positions
        const vertices = [];
        const colors = [];
        const particleColor = new THREE.Color();
        
        for (const flow of flows) {
            const edge = edges[flow.edge];
            const from = centers[edge.from];
            const to = centers[edge.to];
            
            // Interpolate position along edge
            const x = from[0] + (to[0] - from[0]) * flow.position;
            const y = from[1] + (to[1] - from[1]) * flow.position;
            
            vertices.push(x, y, 0);
            
            // Color based on velocity
            const velocityFactor = Math.min(1, Math.abs(flow.velocity) / 2);
            particleColor.setRGB(
                this.flowParticleColor.r * (1 - velocityFactor) + velocityFactor,
                this.flowParticleColor.g * (1 - velocityFactor),
                this.flowParticleColor.b * (1 - velocityFactor)
            );
            
            colors.push(particleColor.r, particleColor.g, particleColor.b);
        }
        
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        
        // Update particle material to use vertex colors
        this.particleMaterial.vertexColors = true;
        
        this.flowParticles = new THREE.Points(geometry, this.particleMaterial);
        this.scene.add(this.flowParticles);
    }
    
    /**
     * Update flow particle positions
     */
    updateFlowParticles(centers, edges, flows) {
        if (!this.flowParticles) return;
        
        const positions = this.flowParticles.geometry.attributes.position.array;
        const colors = this.flowParticles.geometry.attributes.color.array;
        const particleColor = new THREE.Color();
        
        for (let i = 0; i < flows.length; i++) {
            const flow = flows[i];
            const edge = edges[flow.edge];
            const from = centers[edge.from];
            const to = centers[edge.to];
            
            // Update position
            const x = from[0] + (to[0] - from[0]) * flow.position;
            const y = from[1] + (to[1] - from[1]) * flow.position;
            
            positions[i * 3] = x;
            positions[i * 3 + 1] = y;
            
            // Update color based on velocity
            const velocityFactor = Math.min(1, Math.abs(flow.velocity) / 2);
            particleColor.setRGB(
                this.flowParticleColor.r * (1 - velocityFactor) + velocityFactor,
                this.flowParticleColor.g * (1 - velocityFactor),
                this.flowParticleColor.b * (1 - velocityFactor)
            );
            
            colors[i * 3] = particleColor.r;
            colors[i * 3 + 1] = particleColor.g;
            colors[i * 3 + 2] = particleColor.b;
        }
        
        this.flowParticles.geometry.attributes.position.needsUpdate = true;
        this.flowParticles.geometry.attributes.color.needsUpdate = true;
    }
    
    /**
     * Update foam edges visualization when foam structure changes
     */
    updateFoamEdges(centers, edges) {
        if (!this.foamLines) return;
        
        const positions = this.foamLines.geometry.attributes.position.array;
        
        for (let i = 0; i < edges.length; i++) {
            const edge = edges[i];
            const from = centers[edge.from];
            const to = centers[edge.to];
            
            // Update line segment vertices
            positions[i * 6] = from[0];
            positions[i * 6 + 1] = from[1];
            positions[i * 6 + 3] = to[0];
            positions[i * 6 + 4] = to[1];
        }
        
        this.foamLines.geometry.attributes.position.needsUpdate = true;
    }
    
    /**
     * Render the scene
     */
    render() {
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
    
    /**
     * Set visualization colors
     */
    setColors(options) {
        if (options.background) {
            this.backgroundColor.set(options.background);
            this.scene.background = this.backgroundColor;
        }
        
        if (options.triangles) {
            this.triangleColor.set(options.triangles);
            this.triangleMaterial.color = this.triangleColor;
        }
        
        if (options.foamLines) {
            this.foamLineColor.set(options.foamLines);
            this.lineMaterial.color = this.foamLineColor;
        }
        
        if (options.flowParticles) {
            this.flowParticleColor.set(options.flowParticles);
            // Note: this doesn't affect existing particles until updateFlowParticles is called
        }
    }
}

export default FoamRenderer; 