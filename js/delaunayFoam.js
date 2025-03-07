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
        this.foamFlows = [];               // Particles flowing on Voronoi edges
        
        // Center type (circumcenter, barycenter, incenter)
        this.centerType = 'barycenter';
        
        // Dynamics parameters
        this.flowSpeed = 1.0;              // Increased from 0.5 to 1.0
        this.flowStrength = 0.02;          // Increased from 0.01 to 0.02
        this.expansionThreshold = 0.5;
        this.contractionFactor = 0.95;
        this.expansionFactor = 1.05;
        this.equilibriumDistance = 50;
        this.restructureThreshold = 0.05;  // Threshold for re-triangulation
        
        // System stability parameters
        this.maxPointMovement = 10;        // Limit how far a point can move in one step
        this.boundaryPadding = 50;         // Keep points within boundaries
        
        // Performance optimization
        this.bufferFrames = 5;             // Only re-triangulate every N frames
        this.frameCount = 0;
        
        // Initialize points and triangulation
        this.initPoints();
        this.triangulate();
        this.createFoam();
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
     * Initialize flows on Voronoi edges
     */
    initializeFlows() {
        this.foamFlows = [];
        
        // Create multiple flows per edge for better visualization
        for (let i = 0; i < this.voronoiEdges.length; i++) {
            // Add 2-3 particles per edge
            const particlesPerEdge = 2 + Math.floor(Math.random() * 2);
            
            for (let j = 0; j < particlesPerEdge; j++) {
                // Distribute evenly along the edge
                const position = (j + 1) / (particlesPerEdge + 1);
                
                // Ensure good initial velocity
                const velocity = (Math.random() > 0.5 ? 1 : -1) * 
                                  this.flowSpeed * (0.5 + Math.random());
                                  
                this.foamFlows.push({
                    position: position,
                    velocity: velocity,
                    edge: i
                });
            }
        }
        
        console.log(`Created ${this.foamFlows.length} flow particles on ${this.voronoiEdges.length} edges`);
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
     * Updates the simulation by one step
     */
    update(deltaTime) {
        // Update flows along Voronoi edges
        this.updateFlows(deltaTime);
        
        // Transfer flow from Voronoi edges to Delaunay edges
        this.transferVoronoiFlowToDelaunay();
        
        // Update Delaunay points based on edge flows
        const significant = this.updateDelaunayPoints(deltaTime);
        
        // Only re-triangulate if points moved significantly or buffer frames reached
        this.frameCount++;
        if (significant || this.frameCount >= this.bufferFrames) {
            this.triangulate();
            this.createFoam();
            this.frameCount = 0;
        }
    }
    
    /**
     * Update particle flows along Voronoi edges
     */
    updateFlows(deltaTime) {
        for (let i = 0; i < this.foamFlows.length; i++) {
            const flow = this.foamFlows[i];
            const edge = this.voronoiEdges[flow.edge];
            
            // Update position with higher velocity to ensure movement
            flow.position += flow.velocity * deltaTime;
            
            // If reached end of edge, move to next edge
            if (flow.position <= 0 || flow.position >= 1) {
                // Add flow to current edge
                edge.flow += this.flowStrength;
                
                // Determine which vertex we've hit (from or to)
                const vertex = flow.position <= 0 ? edge.from : edge.to;
                
                // Find all Voronoi edges connected to this vertex
                const connectedEdges = [];
                for (let j = 0; j < this.voronoiEdges.length; j++) {
                    if (j !== flow.edge && 
                        (this.voronoiEdges[j].from === vertex || 
                         this.voronoiEdges[j].to === vertex)) {
                        connectedEdges.push(j);
                    }
                }
                
                if (connectedEdges.length > 0) {
                    // Pick a random connected edge
                    const nextEdgeIdx = Math.floor(Math.random() * connectedEdges.length);
                    const nextEdge = connectedEdges[nextEdgeIdx];
                    flow.edge = nextEdge;
                    
                    // Set position at the correct end of the new edge
                    const newEdge = this.voronoiEdges[nextEdge];
                    flow.position = newEdge.from === vertex ? 0 : 1;
                    
                    // Adjust velocity direction based on new position
                    if (flow.position === 0) {
                        flow.velocity = Math.abs(flow.velocity); // Move from 0 toward 1
                    } else {
                        flow.velocity = -Math.abs(flow.velocity); // Move from 1 toward 0
                    }
                    
                    // Randomize velocity a bit to create more interesting flows
                    flow.velocity *= 0.8 + Math.random() * 0.4; // 80% to 120% of original
                } else {
                    // If no connected edges (shouldn't happen in proper mesh), just bounce
                    flow.velocity = -flow.velocity;
                    flow.position = Math.max(0, Math.min(1, flow.position));
                }
            }
        }
    }
    
    /**
     * Transfer flow from Voronoi edges to their corresponding Delaunay edges
     */
    transferVoronoiFlowToDelaunay() {
        for (let i = 0; i < this.voronoiEdges.length; i++) {
            const voronoiEdge = this.voronoiEdges[i];
            
            if (voronoiEdge.delaunayEdge !== undefined && voronoiEdge.flow > 0) {
                // Get the corresponding Delaunay edge
                const delaunayEdge = this.delaunayEdges[voronoiEdge.delaunayEdge];
                
                // Transfer flow with amplification for better visual effect
                const flowTransfer = voronoiEdge.flow * 2.0;
                delaunayEdge.flow += flowTransfer;
                
                // Decay Voronoi edge flow but not completely to maintain motion
                voronoiEdge.flow *= 0.7;
            }
        }
        
        // Amplify flows in connected Delaunay edges for better effect
        for (let i = 0; i < this.delaunayEdges.length; i++) {
            if (this.delaunayEdges[i].flow > 0.05) {
                // Find connected Delaunay edges
                const [p1, p2] = this.delaunayEdges[i].points;
                
                for (let j = 0; j < this.delaunayEdges.length; j++) {
                    if (i === j) continue;
                    
                    const [q1, q2] = this.delaunayEdges[j].points;
                    
                    // Check if edges share a point
                    if (p1 === q1 || p1 === q2 || p2 === q1 || p2 === q2) {
                        // Transfer a small amount of flow to connected edges
                        this.delaunayEdges[j].flow += this.delaunayEdges[i].flow * 0.1;
                    }
                }
            }
        }
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
}

export default DelaunayFoam; 