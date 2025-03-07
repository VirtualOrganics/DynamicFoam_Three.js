import Delaunator from 'delaunator';
import * as THREE from 'three';

/**
 * DelaunayFoam class handles the creation and dynamics of the foam structure
 * based on Delaunay triangulation of points.
 */
class DelaunayFoam {
    constructor(width, height, numPoints = 100) {
        this.width = width;
        this.height = height;
        this.numPoints = numPoints;
        
        // Points for triangulation
        this.points = [];
        this.delaunay = null;
        
        // Delaunay and Voronoi structures
        this.triangleCenters = [];         // Voronoi vertices
        this.delaunayEdges = [];           // Maps [pointIndex1, pointIndex2] => edgeIndex
        this.voronoiEdges = [];            // Stores from/to triangle centers and the corresponding delaunay edge
        this.triangleAdjacency = [];       // For each triangle, lists adjacent triangles
        this.flows = [];               // Particles flowing on Voronoi edges
        
        // Center type (circumcenter, barycenter, incenter)
        this.centerType = 'barycenter';
        
        // Dynamics parameters
        this.flowSpeed = 1.0;              // Speed of particles
        this.flowStrength = 0.02;          // How much particles contract edges
        this.expansionThreshold = 0.5;
        this.contractionFactor = 0.95;
        this.expansionFactor = 1.05;
        this.equilibriumDistance = 50;
        this.restructureThreshold = 0.05;  // Threshold for re-triangulation
        
        // Particle management parameters
        this.particleReleaseInterval = 5.0;  // Time in seconds between particle releases
        this.particleLifetime = 30.0;         // Lifetime of particles in seconds
        this.lastReleaseTime = 0;            // Tracks when we last released particles
        this.maxAngleForTraversal = Math.PI/2; // Maximum angle (90 degrees) for edge traversal
        
        // System stability parameters
        this.maxPointMovement = 10;        // Limit how far a point can move in one step
        this.boundaryPadding = 50;         // Keep points within boundaries
        
        // Performance optimization
        this.bufferFrames = 5;             // Only re-triangulate every N frames
        this.frameCount = 0;
        
        // Time tracking
        this.currentTime = 0;
        this.lastLogTime = 0;
        
        // Initialize points and triangulation
        this.init();
    }
    
    init() {
        // Initialize points
        this.initPoints();
        
        // Perform Delaunay triangulation
        this.triangulate();
        
        // Create Voronoi foam
        this.createFoam();
        
        // Initialize time tracking
        this.currentTime = 0;
        this.lastReleaseTime = 0;
        this.lastLogTime = 0;
    }
    
    /**
     * Generate points using Poisson disc sampling for better distribution
     */
    initPoints() {
        // Simple random distribution for now (will enhance with Poisson later)
        this.points = [];
        for (let i = 0; i < this.numPoints; i++) {
            this.points.push(
                (Math.random() - 0.5) * this.width,
                (Math.random() - 0.5) * this.height
            );
        }
    }
    
    /**
     * Perform Delaunay triangulation on the points and compute edge relationships
     */
    triangulate() {
        // Perform basic Delaunay triangulation
        this.delaunay = new Delaunator(this.points);
        
        // Calculate triangle adjacency - for each triangle, which triangles share an edge with it
        this.calculateTriangleAdjacency();
        
        // Calculate Delaunay edges from triangulation
        this.calculateDelaunayEdges();
    }
    
    /**
     * Calculate adjacency information for each triangle
     */
    calculateTriangleAdjacency() {
        const numTriangles = this.delaunay.triangles.length / 3;
        this.triangleAdjacency = new Array(numTriangles).fill(0).map(() => []);
        
        // Use halfedges to determine adjacency
        for (let i = 0; i < this.delaunay.halfedges.length; i++) {
            const opposite = this.delaunay.halfedges[i];
            if (opposite >= 0) {
                const triangle = Math.floor(i / 3);
                const adjacentTriangle = Math.floor(opposite / 3);
                this.triangleAdjacency[triangle].push(adjacentTriangle);
            }
        }
    }
    
    /**
     * Calculate all Delaunay edges from the triangulation
     */
    calculateDelaunayEdges() {
        const seen = new Set();
        this.delaunayEdges = [];
        
        // Process all triangles to find edges
        for (let t = 0; t < this.delaunay.triangles.length / 3; t++) {
            const i0 = this.delaunay.triangles[3 * t];
            const i1 = this.delaunay.triangles[3 * t + 1];
            const i2 = this.delaunay.triangles[3 * t + 2];
            
            // Add each edge if not already seen
            this.addDelaunayEdgeIfNew(i0, i1, seen);
            this.addDelaunayEdgeIfNew(i1, i2, seen);
            this.addDelaunayEdgeIfNew(i2, i0, seen);
        }
    }
    
    /**
     * Helper to add a Delaunay edge if it hasn't been added yet
     */
    addDelaunayEdgeIfNew(p1, p2, seen) {
        // Ensure consistent ordering (smaller index first)
        const edge = p1 < p2 ? [p1, p2] : [p2, p1];
        const edgeKey = `${edge[0]}-${edge[1]}`;
        
        if (!seen.has(edgeKey)) {
            seen.add(edgeKey);
            // Calculate current length of the edge
            const x1 = this.points[2 * edge[0]];
            const y1 = this.points[2 * edge[0] + 1];
            const x2 = this.points[2 * edge[1]];
            const y2 = this.points[2 * edge[1] + 1];
            const length = Math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2);
            
            this.delaunayEdges.push({
                points: edge,
                length: length,
                flow: 0,
                restLength: length  // Initial rest length is the current length
            });
        }
    }
    
    /**
     * Calculate the circumcenter of a triangle
     */
    calculateCircumcenter(a, b, c) {
        // Get coordinates
        const ax = this.points[2 * a];
        const ay = this.points[2 * a + 1];
        const bx = this.points[2 * b];
        const by = this.points[2 * b + 1];
        const cx = this.points[2 * c];
        const cy = this.points[2 * c + 1];
        
        // Circumcenter calculation
        const D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        if (Math.abs(D) < 1e-10) {
            // Degenerate case, return triangle centroid instead
            return [(ax + bx + cx) / 3, (ay + by + cy) / 3];
        }
        
        const aSq = ax * ax + ay * ay;
        const bSq = bx * bx + by * by;
        const cSq = cx * cx + cy * cy;
        
        const x = (aSq * (by - cy) + bSq * (cy - ay) + cSq * (ay - by)) / D;
        const y = (aSq * (cx - bx) + bSq * (ax - cx) + cSq * (bx - ax)) / D;
        
        return [x, y];
    }
    
    /**
     * Calculate the barycenter (centroid) of a triangle
     */
    calculateBarycenter(a, b, c) {
        const ax = this.points[2 * a];
        const ay = this.points[2 * a + 1];
        const bx = this.points[2 * b];
        const by = this.points[2 * b + 1];
        const cx = this.points[2 * c];
        const cy = this.points[2 * c + 1];
        
        return [(ax + bx + cx) / 3, (ay + by + cy) / 3];
    }
    
    /**
     * Calculate the incenter of a triangle
     */
    calculateIncenter(a, b, c) {
        const ax = this.points[2 * a];
        const ay = this.points[2 * a + 1];
        const bx = this.points[2 * b];
        const by = this.points[2 * b + 1];
        const cx = this.points[2 * c];
        const cy = this.points[2 * c + 1];
        
        // Calculate side lengths
        const a_len = Math.sqrt((bx - cx) ** 2 + (by - cy) ** 2);
        const b_len = Math.sqrt((ax - cx) ** 2 + (ay - cy) ** 2);
        const c_len = Math.sqrt((ax - bx) ** 2 + (ay - by) ** 2);
        
        // Calculate incenter
        const x = (a_len * ax + b_len * bx + c_len * cx) / (a_len + b_len + c_len);
        const y = (a_len * ay + b_len * by + c_len * cy) / (a_len + b_len + c_len);
        
        return [x, y];
    }
    
    /**
     * Create foam structure based on triangle centers
     */
    createFoam() {
        // Clear existing structures
        this.triangleCenters = [];
        this.voronoiEdges = [];
        this.foamFlows = [];
        
        // Calculate centers for each triangle
        for (let i = 0; i < this.delaunay.triangles.length / 3; i++) {
            const a = this.delaunay.triangles[i * 3];
            const b = this.delaunay.triangles[i * 3 + 1];
            const c = this.delaunay.triangles[i * 3 + 2];
            
            // Calculate center based on selected type
            let center;
            if (this.centerType === 'circumcenter') {
                center = this.calculateCircumcenter(a, b, c);
            } else if (this.centerType === 'barycenter') {
                center = this.calculateBarycenter(a, b, c);
            } else {
                center = this.calculateIncenter(a, b, c);
            }
            
            this.triangleCenters.push(center);
        }
        
        // Generate Voronoi edges and map to Delaunay edges
        this.generateVoronoiEdges();
        
        // Initialize flows on edges
        this.initializeFlows();
    }
    
    /**
     * Generate Voronoi edges by connecting adjacent triangle centers and map to Delaunay edges
     */
    generateVoronoiEdges() {
        this.voronoiEdges = [];
        
        // Use triangleAdjacency to create Voronoi edges
        for (let t = 0; t < this.triangleAdjacency.length; t++) {
            for (let adjIdx = 0; adjIdx < this.triangleAdjacency[t].length; adjIdx++) {
                const adjT = this.triangleAdjacency[t][adjIdx];
                
                // Only process each edge once (when t < adjT)
                if (t < adjT) {
                    // Find the shared edge between these triangles
                    const sharedPoints = this.findSharedPoints(t, adjT);
                    
                    if (sharedPoints.length === 2) {
                        // Find corresponding Delaunay edge
                        const [p1, p2] = sharedPoints[0] < sharedPoints[1] ? 
                                        [sharedPoints[0], sharedPoints[1]] : 
                                        [sharedPoints[1], sharedPoints[0]];
                        
                        const delaunayEdgeIdx = this.findDelaunayEdgeIndex(p1, p2);
                        
                        if (delaunayEdgeIdx !== -1) {
                            // Create Voronoi edge
                            this.voronoiEdges.push({
                                from: t,
                                to: adjT,
                                delaunayEdge: delaunayEdgeIdx,
                                flow: 0,
                                length: this.calculateDistance(
                                    this.triangleCenters[t],
                                    this.triangleCenters[adjT]
                                )
                            });
                        }
                    }
                }
            }
        }
    }
    
    /**
     * Find the index of the Delaunay edge between two points
     */
    findDelaunayEdgeIndex(p1, p2) {
        for (let i = 0; i < this.delaunayEdges.length; i++) {
            const edge = this.delaunayEdges[i].points;
            if ((edge[0] === p1 && edge[1] === p2) || 
                (edge[0] === p2 && edge[1] === p1)) {
                return i;
            }
        }
        return -1;
    }
    
    /**
     * Find shared points between two triangles
     */
    findSharedPoints(triangleA, triangleB) {
        const triA = [
            this.delaunay.triangles[triangleA * 3],
            this.delaunay.triangles[triangleA * 3 + 1],
            this.delaunay.triangles[triangleA * 3 + 2]
        ];
        
        const triB = [
            this.delaunay.triangles[triangleB * 3],
            this.delaunay.triangles[triangleB * 3 + 1],
            this.delaunay.triangles[triangleB * 3 + 2]
        ];
        
        return triA.filter(point => triB.includes(point));
    }
    
    /**
     * Calculate distance between two points
     */
    calculateDistance(p1, p2) {
        return Math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2);
    }
    
    /**
     * Initialize flows - now just creates an empty array
     */
    initializeFlows() {
        this.foamFlows = [];
        // We'll create particles during the update cycle instead of initializing them here
    }
    
    /**
     * Update the simulation
     * @param {number} deltaTime - Time since last update in seconds
     */
    update(deltaTime) {
        // Track total simulation time
        this.currentTime += deltaTime;
        
        // Initialize particles if this is the first update and we have edges
        if (this.flows.length === 0 && this.voronoiEdges.length > 0 && this.currentTime > 0.1) {
            console.log("Initializing particles on first update");
            this.releaseParticles(this.currentTime);
            this.lastReleaseTime = this.currentTime;
        }
        
        // Check if it's time to release more particles
        if (this.currentTime - this.lastReleaseTime >= this.particleReleaseInterval) {
            this.releaseParticles(this.currentTime);
            this.lastReleaseTime = this.currentTime;
        }
        
        // Update existing particles
        this.updateParticlePositions(deltaTime);
        
        // Apply dynamics (point movement)
        this.updateDelaunayPoints(deltaTime);
        
        // Re-triangulate if needed
        this.frameCount++;
        if (this.frameCount >= this.bufferFrames) {
            this.frameCount = 0;
            this.triangulate();
            this.createFoam();
        }
    }
    
    /**
     * Update existing particle positions
     */
    updateParticlePositions(deltaTime) {
        // Update existing particles
        for (let i = this.flows.length - 1; i >= 0; i--) {
            const flow = this.flows[i];
            
            // Check if the particle has exceeded its lifetime
            if (this.currentTime - flow.birthTime > this.particleLifetime) {
                // Remove expired particle
                this.flows.splice(i, 1);
                continue;
            }
            
            // Get the edge
            if (flow.edgeIndex < 0 || flow.edgeIndex >= this.voronoiEdges.length) {
                console.warn(`Invalid edge index: ${flow.edgeIndex}`);
                this.flows.splice(i, 1);
                continue;
            }
            
            const edge = this.voronoiEdges[flow.edgeIndex];
            
            // Simple movement along the edge
            flow.t += flow.direction * flow.velocity * deltaTime;
            
            // Check for edge transitions
            if (flow.t <= 0 || flow.t >= 1) {
                // Determine which vertex we've reached
                const vertexIndex = flow.t <= 0 ? edge.from : edge.to;
                
                // Find the next edge to move to
                const nextEdgeIndex = this.findNextEdge(flow.edgeIndex, vertexIndex, flow.direction > 0);
                
                if (nextEdgeIndex === -1) {
                    // No valid next edge, remove the particle
                    this.flows.splice(i, 1);
                    continue;
                }
                
                // Transition to the next edge
                flow.edgeIndex = nextEdgeIndex;
                const newEdge = this.voronoiEdges[nextEdgeIndex];
                
                // Set the t-value for the new edge
                if (newEdge.from === vertexIndex) {
                    flow.t = 0.01; // Start slightly inside the edge
                    flow.direction = 1; // Moving from 'from' to 'to'
                } else {
                    flow.t = 0.99; // Start slightly inside the edge
                    flow.direction = -1; // Moving from 'to' to 'from'
                }
                
                // Apply contraction at the vertex
                this.applyContraction(vertexIndex);
            }
        }
    }
    
    /**
     * Release new particles from edges
     */
    releaseParticles(currentTime) {
        // Ensure we have edges to work with
        if (!this.voronoiEdges || this.voronoiEdges.length === 0) {
            console.warn("No edges available for particles");
            return;
        }
        
        // Release particles on ALL edges initially, then be more selective on subsequent releases
        const firstRelease = this.flows.length === 0;
        const edgesToRelease = firstRelease ? this.voronoiEdges.length : Math.min(50, this.voronoiEdges.length);
        
        // Determine step size based on release type
        const step = firstRelease ? 1 : Math.max(1, Math.floor(this.voronoiEdges.length / edgesToRelease));
        
        console.log(`Releasing particles on ${edgesToRelease} edges (step: ${step})`);
        let count = 0;
        
        for (let i = 0; i < this.voronoiEdges.length; i += step) {
            // Skip edges that have null centers (should not happen but safety check)
            const edge = this.voronoiEdges[i];
            if (!this.triangleCenters[edge.from] || !this.triangleCenters[edge.to]) {
                continue;
            }
            
            // Create two particles per edge, moving in opposite directions
            // First particle (moving forward)
            this.flows.push({
                edgeIndex: i,
                t: 0.5,          // Start in the middle of the edge
                velocity: 0.3,   // Fixed velocity for consistent movement
                direction: 1,    // Forward direction
                birthTime: currentTime
            });
            
            // Second particle (moving backward)
            this.flows.push({
                edgeIndex: i,
                t: 0.5,          // Start in the middle of the edge
                velocity: 0.3,   // Fixed velocity for consistent movement
                direction: -1,   // Backward direction
                birthTime: currentTime
            });
            
            count += 2;
        }
        
        console.log(`Released ${count} particles. Total active: ${this.flows.length}`);
    }
    
    /**
     * Find the next edge to move to from a vertex
     * @param {number} currentEdgeIndex - Index of the current edge
     * @param {number} vertexIndex - Index of the vertex
     * @param {boolean} isForwardDirection - Whether the particle is moving in the forward direction
     * @returns {number} Index of the next edge, or -1 if no valid edge found
     */
    findNextEdge(currentEdgeIndex, vertexIndex, isForwardDirection) {
        // Get all edges connected to this vertex
        const connectedEdges = [];
        
        for (let i = 0; i < this.voronoiEdges.length; i++) {
            if (i === currentEdgeIndex) continue; // Skip current edge
            
            const edge = this.voronoiEdges[i];
            if (edge.from === vertexIndex || edge.to === vertexIndex) {
                connectedEdges.push(i);
            }
        }
        
        if (connectedEdges.length === 0) {
            return -1; // No connected edges
        }
        
        // Get the current edge direction vector
        const currentEdge = this.voronoiEdges[currentEdgeIndex];
        let fromIndex, toIndex;
        
        // Determine the direction vector of the current edge
        if (isForwardDirection) {
            // We were moving from from->to
            fromIndex = currentEdge.from;
            toIndex = currentEdge.to;
        } else {
            // We were moving from to->from
            fromIndex = currentEdge.to;
            toIndex = currentEdge.from;
        }
        
        const fromCenter = this.triangleCenters[fromIndex];
        const toCenter = this.triangleCenters[toIndex];
        
        const currentDirX = toCenter[0] - fromCenter[0];
        const currentDirY = toCenter[1] - fromCenter[1];
        
        // Normalize the direction vector
        const currentLength = Math.sqrt(currentDirX * currentDirX + currentDirY * currentDirY);
        const currentNormX = currentDirX / currentLength;
        const currentNormY = currentDirY / currentLength;
        
        // Find the edge with the smallest angle to the current edge
        let bestEdgeIndex = -1;
        let smallestAngle = Math.PI; // Initialize to 180 degrees
        
        for (const nextEdgeIndex of connectedEdges) {
            const nextEdge = this.voronoiEdges[nextEdgeIndex];
            
            // Determine the other vertex of the next edge
            const otherVertexIndex = (nextEdge.from === vertexIndex) ? nextEdge.to : nextEdge.from;
            
            // Get the direction vector of the next edge
            const nextCenter = this.triangleCenters[otherVertexIndex];
            const vertexCenter = this.triangleCenters[vertexIndex];
            
            const nextDirX = nextCenter[0] - vertexCenter[0];
            const nextDirY = nextCenter[1] - vertexCenter[1];
            
            // Normalize the next direction vector
            const nextLength = Math.sqrt(nextDirX * nextDirX + nextDirY * nextDirY);
            const nextNormX = nextDirX / nextLength;
            const nextNormY = nextDirY / nextLength;
            
            // Calculate the angle between the vectors (dot product)
            const dotProduct = currentNormX * nextNormX + currentNormY * nextNormY;
            
            // Convert to angle (in radians)
            const angle = Math.acos(Math.min(1, Math.max(-1, dotProduct)));
            
            // Only consider edges that form less than 90 degrees
            if (angle < Math.PI / 2 && angle < smallestAngle) {
                smallestAngle = angle;
                bestEdgeIndex = nextEdgeIndex;
            }
        }
        
        // If no edge is less than 90 degrees, just pick the one with the smallest angle
        if (bestEdgeIndex === -1 && connectedEdges.length > 0) {
            bestEdgeIndex = connectedEdges[0];
        }
        
        return bestEdgeIndex;
    }
    
    /**
     * Calculate the angle between three consecutive points
     */
    calculateAngle(p1, p2, p3) {
        const a1 = Math.atan2(p1[1] - p2[1], p1[0] - p2[0]);
        const a2 = Math.atan2(p3[1] - p2[1], p3[0] - p2[0]);
        
        let angle = a2 - a1;
        while (angle < 0) angle += 2 * Math.PI;
        while (angle > 2 * Math.PI) angle -= 2 * Math.PI;
        
        return Math.min(angle, 2 * Math.PI - angle);
    }
    
    /**
     * Update Delaunay points based on edge flows
     * Returns true if significant changes occurred
     */
    updateDelaunayPoints(deltaTime) {
        // Calculate point movements based on edge flows
        const movements = new Array(this.points.length / 2).fill(0).map(() => [0, 0]);
        let maxMovement = 0;
        
        // For each Delaunay edge with flow
        for (let i = 0; i < this.delaunayEdges.length; i++) {
            const edge = this.delaunayEdges[i];
            
            // Skip edges with little flow
            if (Math.abs(edge.flow) < 0.01) {
                // Decay flow over time for all edges
                edge.flow *= 0.99;
                continue;
            }
            
            const [p1, p2] = edge.points;
            
            // Calculate current distance
            const p1x = this.points[2 * p1];
            const p1y = this.points[2 * p1 + 1];
            const p2x = this.points[2 * p2];
            const p2y = this.points[2 * p2 + 1];
            
            const currentDistance = Math.sqrt((p2x - p1x) ** 2 + (p2y - p1y) ** 2);
            
            // Target distance based on flow
            let targetDistance = edge.restLength; // Use rest length as base
            
            if (edge.flow > this.expansionThreshold) {
                // Expansion for high flow
                targetDistance *= this.expansionFactor;
            } else {
                // Contraction for normal flow
                targetDistance *= this.contractionFactor;
            }
            
            // Calculate force
            const factor = (targetDistance / currentDistance - 1) * 0.1;
            const dx = (p2x - p1x) * factor;
            const dy = (p2y - p1y) * factor;
            
            // Accumulate movements
            movements[p1][0] -= dx;
            movements[p1][1] -= dy;
            movements[p2][0] += dx;
            movements[p2][1] += dy;
            
            // Track max movement for re-triangulation decision
            const movementMagnitude = Math.sqrt(dx * dx + dy * dy);
            maxMovement = Math.max(maxMovement, movementMagnitude);
            
            // Decay flow over time
            edge.flow *= 0.95;
        }
        
        // Apply movements to points with limits
        for (let i = 0; i < movements.length; i++) {
            // Limit maximum movement
            const movementMagnitude = Math.sqrt(movements[i][0] ** 2 + movements[i][1] ** 2);
            if (movementMagnitude > this.maxPointMovement) {
                const scale = this.maxPointMovement / movementMagnitude;
                movements[i][0] *= scale;
                movements[i][1] *= scale;
            }
            
            // Apply movement
            this.points[2 * i] += movements[i][0];
            this.points[2 * i + 1] += movements[i][1];
            
            // Keep points within boundaries
            this.points[2 * i] = Math.max(-this.width/2 + this.boundaryPadding, 
                              Math.min(this.width/2 - this.boundaryPadding, this.points[2 * i]));
            this.points[2 * i + 1] = Math.max(-this.height/2 + this.boundaryPadding, 
                                Math.min(this.height/2 - this.boundaryPadding, this.points[2 * i + 1]));
        }
        
        // Return whether significant movement occurred
        return maxMovement > this.restructureThreshold;
    }
    
    /**
     * Return the geometry data for Three.js visualization
     */
    getGeometryData() {
        return {
            points: this.points,
            triangles: Array.from(this.delaunay.triangles),
            centers: this.triangleCenters,
            edges: this.voronoiEdges,
            flows: this.foamFlows,
            delaunayEdges: this.delaunayEdges
        };
    }
    
    /**
     * Change the center type and update the foam
     */
    setCenterType(type) {
        if (['circumcenter', 'barycenter', 'incenter'].includes(type)) {
            this.centerType = type;
            this.createFoam();
            return true;
        }
        return false;
    }
    
    /**
     * Set the flow dynamics parameters
     */
    setDynamicsParams(params) {
        if (params.flowSpeed !== undefined) this.flowSpeed = params.flowSpeed;
        if (params.flowStrength !== undefined) this.flowStrength = params.flowStrength;
        if (params.expansionThreshold !== undefined) this.expansionThreshold = params.expansionThreshold;
        if (params.contractionFactor !== undefined) this.contractionFactor = params.contractionFactor;
        if (params.expansionFactor !== undefined) this.expansionFactor = params.expansionFactor;
        if (params.equilibriumDistance !== undefined) this.equilibriumDistance = params.equilibriumDistance;
        if (params.restructureThreshold !== undefined) this.restructureThreshold = params.restructureThreshold;
        if (params.maxAngleForTraversal !== undefined) this.maxAngleForTraversal = params.maxAngleForTraversal * Math.PI / 180;
        if (params.particleReleaseInterval !== undefined) this.particleReleaseInterval = params.particleReleaseInterval;
        if (params.particleLifetime !== undefined) this.particleLifetime = params.particleLifetime;
        if (params.maxPointMovement !== undefined) this.maxPointMovement = params.maxPointMovement;
        if (params.boundaryPadding !== undefined) this.boundaryPadding = params.boundaryPadding;
        if (params.bufferFrames !== undefined) this.bufferFrames = params.bufferFrames;
        if (params.currentTime !== undefined) this.currentTime = params.currentTime;
        if (params.lastReleaseTime !== undefined) this.lastReleaseTime = params.lastReleaseTime;
        if (params.lastLogTime !== undefined) this.lastLogTime = params.lastLogTime;
        this.createFoam();
    }
    
    /**
     * Reset the simulation with new points
     */
    reset(numPoints = null) {
        if (numPoints !== null) this.numPoints = numPoints;
        this.initPoints();
        this.triangulate();
        this.createFoam();
    }
    
    /**
     * Apply contraction effect at a vertex
     * @param {number} vertexIndex - Index of the vertex
     */
    applyContraction(vertexIndex) {
        if (vertexIndex < 0 || vertexIndex >= this.triangleCenters.length) {
            return;
        }
        
        // Contract Delaunay edges connected to this vertex
        for (let i = 0; i < this.delaunayEdges.length; i++) {
            const edge = this.delaunayEdges[i];
            if (edge.from === vertexIndex || edge.to === vertexIndex) {
                // Apply contraction effect
                edge.flow += this.flowStrength * 0.5;
            }
        }
    }
}

export default DelaunayFoam; 