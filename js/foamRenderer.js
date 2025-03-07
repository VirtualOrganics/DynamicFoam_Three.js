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
        this.triangles = null;         // Delaunay triangulation
        this.foamLines = null;         // Voronoi edges
        this.flowParticles = null;     // Flow particles
        this.delaunayLines = null;     // Delaunay edges
        this.acuteCorners = null;      // Markers for acute angles
        
        // Materials
        this.triangleMaterial = null;
        this.lineMaterial = null;
        this.particleMaterial = null;
        this.delaunayMaterial = null;
        
        // Colors
        this.backgroundColor = new THREE.Color(0x111111);
        this.triangleColor = new THREE.Color(0x444444);
        this.foamLineColor = new THREE.Color(0x2288ff);
        this.flowParticleColor = new THREE.Color(0xff8822);
        this.delaunayLineColor = new THREE.Color(0x7a7a7a);
        this.acuteCornerColor = new THREE.Color(0xff0000);
        
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
            size: 6,          // Increased from 3 to 6
            sizeAttenuation: true,
            transparent: true,
            opacity: 0.8
        });
        
        this.delaunayMaterial = new THREE.LineBasicMaterial({
            color: this.delaunayLineColor,
            linewidth: 1,
            transparent: true,
            opacity: 0.7
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
     * Update foam edges (Voronoi) visualization when foam structure changes
     */
    updateFoamEdges(centers, edges) {
        if (!this.foamLines) return;
        
        const positions = this.foamLines.geometry.attributes.position.array;
        const colors = this.foamLines.geometry.attributes.color.array;
        
        for (let i = 0; i < edges.length; i++) {
            const edge = edges[i];
            const from = centers[edge.from];
            const to = centers[edge.to];
            
            // Update line segment vertices
            positions[i * 6] = from[0];
            positions[i * 6 + 1] = from[1];
            positions[i * 6 + 3] = to[0];
            positions[i * 6 + 4] = to[1];
            
            // Update color based on flow
            const flowIntensity = Math.min(1, Math.abs(edge.flow) * 10);
            const color = new THREE.Color(this.foamLineColor);
            color.lerp(new THREE.Color(0x00ffff), flowIntensity);
            
            colors[i * 6] = color.r;
            colors[i * 6 + 1] = color.g;
            colors[i * 6 + 2] = color.b;
            colors[i * 6 + 3] = color.r;
            colors[i * 6 + 4] = color.g;
            colors[i * 6 + 5] = color.b;
        }
        
        this.foamLines.geometry.attributes.position.needsUpdate = true;
        this.foamLines.geometry.attributes.color.needsUpdate = true;
    }
    
    /**
     * Update Delaunay edges visualization
     */
    updateDelaunayEdges(points, delaunayEdges) {
        if (!this.delaunayLines) return;
        
        const positions = this.delaunayLines.geometry.attributes.position.array;
        const colors = this.delaunayLines.geometry.attributes.color.array;
        
        for (let i = 0; i < delaunayEdges.length; i++) {
            const edge = delaunayEdges[i];
            const p1 = edge.points[0];
            const p2 = edge.points[1];
            
            // Update line segment vertices
            positions[i * 6] = points[2 * p1];
            positions[i * 6 + 1] = points[2 * p1 + 1];
            positions[i * 6 + 3] = points[2 * p2];
            positions[i * 6 + 4] = points[2 * p2 + 1];
            
            // Update color based on flow
            const flowIntensity = Math.min(1, Math.abs(edge.flow) * 10);
            const color = new THREE.Color(this.delaunayLineColor);
            color.lerp(new THREE.Color(0xff0000), flowIntensity);
            
            colors[i * 6] = color.r;
            colors[i * 6 + 1] = color.g;
            colors[i * 6 + 2] = color.b;
            colors[i * 6 + 3] = color.r;
            colors[i * 6 + 4] = color.g;
            colors[i * 6 + 5] = color.b;
        }
        
        this.delaunayLines.geometry.attributes.position.needsUpdate = true;
        this.delaunayLines.geometry.attributes.color.needsUpdate = true;
    }
    
    /**
     * Update FoamData - override to check particle rendering
     */
    updateFoamData(foamData) {
        // Store current particle visibility state
        const particlesVisible = this.flowParticles ? this.flowParticles.visible : true;
        
        // Print debugging info about incoming data
        console.log("updateFoamData called with:", {
            pointsLength: foamData.points ? foamData.points.length : 0,
            trianglesLength: foamData.triangles ? foamData.triangles.length : 0,
            delaunayEdgesLength: foamData.delaunayEdges ? foamData.delaunayEdges.length : 0,
            centersLength: foamData.centers ? foamData.centers.length : 0,
            edgesLength: foamData.edges ? foamData.edges.length : 0,
            flowsLength: foamData.flows ? foamData.flows.length : 0
        });
        
        if (foamData.flows && foamData.flows.length > 0) {
            // Log a sample flow for debugging
            console.log("Sample flow:", JSON.stringify(foamData.flows[0]));
        }
        
        // Remove previous objects except particles
        if (this.triangles) {
            this.scene.remove(this.triangles);
            this.triangles = null;
        }
        if (this.foamLines) {
            this.scene.remove(this.foamLines);
            this.foamLines = null;
        }
        if (this.delaunayLines) {
            this.scene.remove(this.delaunayLines);
            this.delaunayLines = null;
        }
        if (this.acuteCorners) {
            this.scene.remove(this.acuteCorners);
            this.acuteCorners = null;
        }
        
        // Create delaunay edges visualization
        this.createDelaunayEdges(foamData.points, foamData.delaunayEdges);
        
        // Create foam edges (Voronoi)
        this.createFoamEdges(foamData.centers, foamData.edges);
        
        // Create acute angle markers
        this.createAcuteCorners(foamData.centers, foamData.edges);
        
        // Update flow particles (don't recreate unless necessary)
        this.updateFlowParticles(foamData.centers, foamData.edges, foamData.flows);
        
        // Restore particle visibility
        if (this.flowParticles) {
            this.flowParticles.visible = particlesVisible;
            console.log(`Particle visibility set to ${particlesVisible}, particle count: ${this.flowParticles.children ? this.flowParticles.children.length : 0}`);
        } else {
            console.warn("flowParticles is null after updateFlowParticles call");
        }
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
     * Create visualization of Delaunay edges
     */
    createDelaunayEdges(points, delaunayEdges) {
        const geometry = new THREE.BufferGeometry();
        
        // Create line segments from Delaunay edges
        const vertices = [];
        for (const edge of delaunayEdges) {
            const p1 = edge.points[0];
            const p2 = edge.points[1];
            
            vertices.push(
                points[2 * p1], points[2 * p1 + 1], 0,
                points[2 * p2], points[2 * p2 + 1], 0
            );
        }
        
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        
        // Create line colors based on flow
        const colors = [];
        for (const edge of delaunayEdges) {
            // Color intensity based on flow
            const flowIntensity = Math.min(1, Math.abs(edge.flow) * 10);
            const color = new THREE.Color(this.delaunayLineColor);
            color.lerp(new THREE.Color(0xff0000), flowIntensity);
            
            colors.push(color.r, color.g, color.b);
            colors.push(color.r, color.g, color.b);
        }
        
        geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        
        const material = new THREE.LineBasicMaterial({
            vertexColors: true,
            linewidth: 1,
            transparent: true,
            opacity: 0.7
        });
        
        this.delaunayLines = new THREE.LineSegments(geometry, material);
        this.scene.add(this.delaunayLines);
    }
    
    /**
     * Create visualization of foam edges (Voronoi edges)
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
        
        // Create line colors based on flow
        const colors = [];
        for (const edge of edges) {
            // Color intensity based on flow
            const flowIntensity = Math.min(1, Math.abs(edge.flow) * 10);
            const color = new THREE.Color(this.foamLineColor);
            color.lerp(new THREE.Color(0x00ffff), flowIntensity);
            
            colors.push(color.r, color.g, color.b);
            colors.push(color.r, color.g, color.b);
        }
        
        geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        
        const material = new THREE.LineBasicMaterial({
            vertexColors: true,
            linewidth: 1
        });
        
        this.foamLines = new THREE.LineSegments(geometry, material);
        this.scene.add(this.foamLines);
    }
    
    /**
     * Create flow particles visualization using spheres
     */
    createFlowParticles(centers, edges, flows) {
        // Remove old particles if they exist
        if (this.flowParticles) {
            // If it's a group, remove all children
            if (this.flowParticles.isGroup) {
                while (this.flowParticles.children.length > 0) {
                    const mesh = this.flowParticles.children[0];
                    mesh.geometry.dispose();
                    mesh.material.dispose();
                    this.flowParticles.remove(mesh);
                }
            }
            
            this.scene.remove(this.flowParticles);
            this.flowParticles = null;
        }
        
        if (!flows || flows.length === 0) {
            console.log("No flow particles to visualize");
            return;
        }
        
        console.log(`Creating flow particles for ${flows.length} particles`);
        
        // Create a new group to hold all particle meshes
        this.flowParticles = new THREE.Group();
        
        // Create sphere geometry and material once for reuse
        const sphereGeometry = new THREE.SphereGeometry(6, 8, 8); // Smaller radius
        const sphereMaterial = new THREE.MeshBasicMaterial({
            color: 0xffffff,  // Pure white
            transparent: false,
            opacity: 1.0
        });
        
        let validParticleCount = 0;
        
        // Create a mesh for each particle
        for (let i = 0; i < flows.length; i++) {
            const flow = flows[i];
            
            if (!flow || flow.edgeIndex === undefined || flow.edgeIndex < 0 || flow.edgeIndex >= edges.length) {
                continue;
            }
            
            const edge = edges[flow.edgeIndex];
            if (!edge || edge.from === undefined || edge.to === undefined) {
                continue;
            }
            
            const fromCenter = centers[edge.from];
            const toCenter = centers[edge.to];
            
            if (!fromCenter || !toCenter) {
                continue;
            }
            
            // Calculate position based on t parameter
            const t = flow.t !== undefined ? flow.t : 0.5;
            const x = fromCenter[0] + t * (toCenter[0] - fromCenter[0]);
            const y = fromCenter[1] + t * (toCenter[1] - fromCenter[1]);
            
            // Create a sphere mesh for this particle
            const particleMesh = new THREE.Mesh(sphereGeometry, sphereMaterial.clone());
            particleMesh.position.set(x, y, 5); // Higher z-offset
            particleMesh.userData.flowIndex = i; // Store original index
            
            // Add to the group
            this.flowParticles.add(particleMesh);
            validParticleCount++;
        }
        
        if (validParticleCount === 0) {
            console.warn("No valid particles to render");
            return;
        }
        
        // Add the group to the scene
        this.scene.add(this.flowParticles);
        console.log(`Added ${validParticleCount} flow particle meshes to scene`);
    }
    
    /**
     * Update flow particles positions based on flows
     */
    updateFlowParticles(centers, edges, flows) {
        // If particles don't exist or flows array size has drastically changed, recreate them
        const shouldRecreate = !this.flowParticles || 
                               !flows || 
                               !this.flowParticles.isGroup ||
                               (flows.length > 0 && Math.abs(this.flowParticles.children.length - flows.length) > 20);
        
        if (shouldRecreate) {
            this.createFlowParticles(centers, edges, flows);
            return;
        }
        
        if (flows.length === 0) {
            return; // No particles to update
        }
        
        // Update positions of existing meshes, add or remove as needed
        const minCount = Math.min(flows.length, this.flowParticles.children.length);
        
        // Update existing particles
        for (let i = 0; i < minCount; i++) {
            const flow = flows[i];
            const mesh = this.flowParticles.children[i];
            
            if (!flow || flow.edgeIndex === undefined || flow.edgeIndex < 0 || flow.edgeIndex >= edges.length) {
                mesh.visible = false;
                continue;
            }
            
            const edge = edges[flow.edgeIndex];
            if (!edge || edge.from === undefined || edge.to === undefined) {
                mesh.visible = false;
                continue;
            }
            
            const fromCenter = centers[edge.from];
            const toCenter = centers[edge.to];
            
            if (!fromCenter || !toCenter) {
                mesh.visible = false;
                continue;
            }
            
            // Calculate new position
            const t = flow.t !== undefined ? flow.t : 0.5;
            mesh.position.x = fromCenter[0] + t * (toCenter[0] - fromCenter[0]);
            mesh.position.y = fromCenter[1] + t * (toCenter[1] - fromCenter[1]);
            mesh.position.z = 5; // Keep above the plane
            mesh.visible = true;
        }
        
        // Handle case where we have more flows than meshes
        if (flows.length > this.flowParticles.children.length) {
            const sphereGeometry = new THREE.SphereGeometry(6, 8, 8);
            const sphereMaterial = new THREE.MeshBasicMaterial({
                color: 0xffffff,
                transparent: false,
                opacity: 1.0
            });
            
            // Create additional particles
            for (let i = this.flowParticles.children.length; i < flows.length; i++) {
                const flow = flows[i];
                
                if (!flow || flow.edgeIndex === undefined || flow.edgeIndex < 0 || flow.edgeIndex >= edges.length) {
                    continue;
                }
                
                const edge = edges[flow.edgeIndex];
                if (!edge || !centers[edge.from] || !centers[edge.to]) {
                    continue;
                }
                
                const fromCenter = centers[edge.from];
                const toCenter = centers[edge.to];
                
                // Calculate position
                const t = flow.t !== undefined ? flow.t : 0.5;
                const x = fromCenter[0] + t * (toCenter[0] - fromCenter[0]);
                const y = fromCenter[1] + t * (toCenter[1] - fromCenter[1]);
                
                // Create new mesh
                const particleMesh = new THREE.Mesh(sphereGeometry, sphereMaterial.clone());
                particleMesh.position.set(x, y, 5);
                
                // Add to the group
                this.flowParticles.add(particleMesh);
            }
        }
        // Handle case where we have more meshes than flows
        else if (flows.length < this.flowParticles.children.length) {
            // Hide excess meshes
            for (let i = flows.length; i < this.flowParticles.children.length; i++) {
                this.flowParticles.children[i].visible = false;
            }
        }
    }
    
    /**
     * Create markers for acute angle corners
     */
    createAcuteCorners(centers, edges) {
        // Remove old markers if they exist
        if (this.acuteCorners) {
            this.scene.remove(this.acuteCorners);
            this.acuteCorners = null;
        }
        
        // Find vertices with acute angles
        const acuteVertices = this.findAcuteAngleVertices(centers, edges);
        
        if (acuteVertices.length === 0) {
            console.log("No acute angle vertices found");
            return;
        }
        
        console.log(`Found ${acuteVertices.length} acute angle vertices`);
        
        // Create a geometry for the markers
        const geometry = new THREE.BufferGeometry();
        const vertices = [];
        
        for (const vertex of acuteVertices) {
            vertices.push(vertex.x, vertex.y, 2); // Same z-offset as particles
        }
        
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        
        // Create a material for the markers (red)
        const material = new THREE.PointsMaterial({
            color: this.acuteCornerColor,
            size: 15,
            sizeAttenuation: true,
            transparent: true,
            opacity: 0.8,
            map: this.createCornerTexture()
        });
        
        this.acuteCorners = new THREE.Points(geometry, material);
        this.scene.add(this.acuteCorners);
    }
    
    /**
     * Create a different texture for corner markers (X shape)
     */
    createCornerTexture() {
        const canvas = document.createElement('canvas');
        canvas.width = 64;
        canvas.height = 64;
        
        const context = canvas.getContext('2d');
        const centerX = canvas.width / 2;
        const centerY = canvas.height / 2;
        const size = canvas.width / 2.5;
        
        // Clear canvas
        context.clearRect(0, 0, canvas.width, canvas.height);
        
        // Draw X shape
        context.strokeStyle = 'white';
        context.lineWidth = 4;
        
        context.beginPath();
        context.moveTo(centerX - size, centerY - size);
        context.lineTo(centerX + size, centerY + size);
        context.stroke();
        
        context.beginPath();
        context.moveTo(centerX + size, centerY - size);
        context.lineTo(centerX - size, centerY + size);
        context.stroke();
        
        // Create radial gradient for glow
        const gradient = context.createRadialGradient(
            centerX, centerY, 0,
            centerX, centerY, size
        );
        gradient.addColorStop(0, 'rgba(255, 0, 0, 0.8)');
        gradient.addColorStop(1, 'rgba(255, 0, 0, 0)');
        
        // Draw glow
        context.fillStyle = gradient;
        context.beginPath();
        context.arc(centerX, centerY, size, 0, Math.PI * 2, false);
        context.fill();
        
        // Create texture
        const texture = new THREE.CanvasTexture(canvas);
        texture.needsUpdate = true;
        return texture;
    }
    
    /**
     * Find vertices where adjacent edges form acute angles
     */
    findAcuteAngleVertices(centers, edges) {
        const acuteVertices = [];
        const vertexMap = new Map(); // Map vertices to their connected edges
        
        // Map edges to vertices
        for (let i = 0; i < edges.length; i++) {
            const edge = edges[i];
            
            // Add edge to 'from' vertex
            if (!vertexMap.has(edge.from)) {
                vertexMap.set(edge.from, []);
            }
            vertexMap.get(edge.from).push(i);
            
            // Add edge to 'to' vertex
            if (!vertexMap.has(edge.to)) {
                vertexMap.set(edge.to, []);
            }
            vertexMap.get(edge.to).push(i);
        }
        
        // Check all vertices for acute angles
        for (const [vertexIdx, connectedEdges] of vertexMap.entries()) {
            // Need at least 2 edges to form an angle
            if (connectedEdges.length < 2) continue;
            
            const vertexPos = centers[vertexIdx];
            let hasAcuteAngle = false;
            
            // Check all pairs of edges
            for (let i = 0; i < connectedEdges.length; i++) {
                const edge1 = edges[connectedEdges[i]];
                const edge1Endpoint = edge1.from === vertexIdx ? edge1.to : edge1.from;
                const edge1Dir = [
                    centers[edge1Endpoint][0] - vertexPos[0],
                    centers[edge1Endpoint][1] - vertexPos[1]
                ];
                
                for (let j = i + 1; j < connectedEdges.length; j++) {
                    const edge2 = edges[connectedEdges[j]];
                    const edge2Endpoint = edge2.from === vertexIdx ? edge2.to : edge2.from;
                    const edge2Dir = [
                        centers[edge2Endpoint][0] - vertexPos[0],
                        centers[edge2Endpoint][1] - vertexPos[1]
                    ];
                    
                    // Normalize vectors
                    const len1 = Math.sqrt(edge1Dir[0] * edge1Dir[0] + edge1Dir[1] * edge1Dir[1]);
                    const len2 = Math.sqrt(edge2Dir[0] * edge2Dir[0] + edge2Dir[1] * edge2Dir[1]);
                    
                    edge1Dir[0] /= len1;
                    edge1Dir[1] /= len1;
                    edge2Dir[0] /= len2;
                    edge2Dir[1] /= len2;
                    
                    // Calculate dot product
                    const dotProduct = edge1Dir[0] * edge2Dir[0] + edge1Dir[1] * edge2Dir[1];
                    
                    // Acute angle if dot product > 0 (angle < 90 degrees)
                    if (dotProduct > 0) {
                        hasAcuteAngle = true;
                        break;
                    }
                }
                
                if (hasAcuteAngle) break;
            }
            
            if (hasAcuteAngle) {
                acuteVertices.push({ x: vertexPos[0], y: vertexPos[1] });
            }
        }
        
        return acuteVertices;
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
        
        if (options.delaunayLines) {
            this.delaunayLineColor.set(options.delaunayLines);
            this.delaunayMaterial.color = this.delaunayLineColor;
        }
    }
    
    /**
     * Toggle visibility of different components
     */
    setVisibility(options) {
        if (this.triangles && options.showTriangulation !== undefined) {
            this.triangles.visible = options.showTriangulation;
        }
        
        if (this.foamLines && options.showFoamEdges !== undefined) {
            this.foamLines.visible = options.showFoamEdges;
        }
        
        if (this.flowParticles && options.showFlowParticles !== undefined) {
            this.flowParticles.visible = options.showFlowParticles;
        }
        
        if (this.delaunayLines && options.showDelaunayEdges !== undefined) {
            this.delaunayLines.visible = options.showDelaunayEdges;
        }
        
        if (this.acuteCorners && options.showAcuteCorners !== undefined) {
            this.acuteCorners.visible = options.showAcuteCorners;
        }
    }
    
    /**
     * Toggle dark/light mode
     */
    setDarkMode(isDarkMode) {
        if (isDarkMode) {
            this.backgroundColor = new THREE.Color(0x111111);
            this.triangleColor = new THREE.Color(0x444444);
            this.foamLineColor = new THREE.Color(0x2288ff);
            this.delaunayLineColor = new THREE.Color(0x7a7a7a);
        } else {
            this.backgroundColor = new THREE.Color(0xf0f0f0);
            this.triangleColor = new THREE.Color(0xdddddd);
            this.foamLineColor = new THREE.Color(0x0066cc);
            this.delaunayLineColor = new THREE.Color(0x999999);
        }
        
        // Update the background color
        this.renderer.setClearColor(this.backgroundColor);
        
        // Update materials if they exist
        if (this.triangleMaterial) {
            this.triangleMaterial.color = this.triangleColor;
        }
        
        if (this.lineMaterial) {
            this.lineMaterial.color = this.foamLineColor;
        }
        
        if (this.delaunayMaterial) {
            this.delaunayMaterial.color = this.delaunayLineColor;
        }
    }
}

export default FoamRenderer; 