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
        
        // Foam structures
        this.triangleCenters = [];
        this.foamEdges = [];
        this.foamFlows = [];
        
        // Center type (circumcenter, barycenter, incenter)
        this.centerType = 'circumcenter';
        
        // Dynamics parameters
        this.flowSpeed = 0.5;
        this.flowStrength = 0.01;
        this.expansionThreshold = 0.5;
        this.contractionFactor = 0.95;
        this.expansionFactor = 1.05;
        this.equilibriumDistance = 50;
        
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
     * Perform Delaunay triangulation on the points
     */
    triangulate() {
        this.delaunay = new Delaunator(this.points);
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
        this.triangleCenters = [];
        
        // Calculate centers for each triangle
        for (let i = 0; i < this.delaunay.triangles.length; i += 3) {
            const a = this.delaunay.triangles[i];
            const b = this.delaunay.triangles[i + 1];
            const c = this.delaunay.triangles[i + 2];
            
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
        
        // Generate foam edges by connecting adjacent triangle centers
        this.generateFoamEdges();
        
        // Initialize flows on edges
        this.initializeFlows();
    }
    
    /**
     * Generate foam edges by connecting adjacent triangle centers
     */
    generateFoamEdges() {
        this.foamEdges = [];
        const halfedges = this.delaunay.halfedges;
        
        for (let i = 0; i < this.delaunay.triangles.length / 3; i++) {
            for (let j = 0; j < 3; j++) {
                const edgeIndex = 3 * i + j;
                const oppositeEdge = halfedges[edgeIndex];
                
                // Skip if we've already processed this edge or it's an external edge
                if (oppositeEdge === -1 || oppositeEdge < edgeIndex) continue;
                
                const triangleA = i;
                const triangleB = Math.floor(oppositeEdge / 3);
                
                // Add the edge between the centers of these triangles
                this.foamEdges.push({
                    from: triangleA,
                    to: triangleB,
                    flow: 0,
                    length: this.calculateDistance(
                        this.triangleCenters[triangleA],
                        this.triangleCenters[triangleB]
                    )
                });
            }
        }
    }
    
    /**
     * Calculate distance between two points
     */
    calculateDistance(p1, p2) {
        return Math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2);
    }
    
    /**
     * Initialize flows on all edges to zero
     */
    initializeFlows() {
        this.foamFlows = [];
        for (let i = 0; i < this.foamEdges.length; i++) {
            this.foamFlows.push({
                position: 0.5, // Start in middle of edge
                velocity: this.flowSpeed * (Math.random() > 0.5 ? 1 : -1), // Use consistent speed with random direction
                edge: i
            });
        }
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
     * Update particle flows along edges
     */
    updateFlows(deltaTime) {
        for (let i = 0; i < this.foamFlows.length; i++) {
            const flow = this.foamFlows[i];
            const edge = this.foamEdges[flow.edge];
            
            // Update velocity based on flowSpeed parameter (keep direction, update magnitude)
            flow.velocity = Math.sign(flow.velocity) * this.flowSpeed;
            
            // Previous position (to detect crossing a Delaunay edge)
            const prevPosition = flow.position;
            
            // Update position - this is what makes particles move along edges
            flow.position += flow.velocity * deltaTime;
            
            // Detect if particle crosses the middle of the Voronoi edge (where Delaunay edge crosses)
            const crossedMiddle = (prevPosition < 0.5 && flow.position >= 0.5) || 
                                 (prevPosition >= 0.5 && flow.position < 0.5);
            
            if (crossedMiddle) {
                // When crossing the middle of a Voronoi edge, apply flow to corresponding Delaunay edge
                const triangleA = edge.from;
                const triangleB = edge.to;
                
                // Find the shared edge between these triangles (the Delaunay edge)
                const sharedPoints = this.findSharedPoints(triangleA, triangleB);
                
                if (sharedPoints.length == 2) {
                    // Apply flow to this Delaunay edge
                    const p1 = sharedPoints[0];
                    const p2 = sharedPoints[1];
                    
                    // Apply flow to this edge - this is what affects the Delaunay mesh
                    this.addFlowToEdge(p1, p2, this.flowStrength);
                }
            }
            
            // If particle reaches end of edge, transition to next edge
            if (flow.position >= 1.0) {
                flow.position = 1.0; // Clamp to edge end
                
                // Find the next edge connected to the "to" vertex
                const connectedEdges = this.foamEdges.filter((e, idx) => 
                    (e.from === edge.to || e.to === edge.to) && 
                    idx !== flow.edge
                );
                
                if (connectedEdges.length > 0) {
                    // Pick a connected edge to move to
                    const nextEdgeIdx = Math.floor(Math.random() * connectedEdges.length);
                    const nextEdge = connectedEdges[nextEdgeIdx];
                    flow.edge = this.foamEdges.indexOf(nextEdge);
                    
                    // Start at the correct end of the new edge
                    flow.position = nextEdge.from === edge.to ? 0.0 : 1.0;
                    
                    // Make sure the velocity is in the correct direction to move away from junction
                    if ((flow.position === 0.0 && flow.velocity < 0) || 
                        (flow.position === 1.0 && flow.velocity > 0)) {
                        flow.velocity = -flow.velocity;
                    }
                } else {
                    // No connected edges, just reverse direction
                    flow.velocity = -flow.velocity;
                    flow.position = 0.99; // Move slightly back to avoid getting stuck
                }
            } else if (flow.position <= 0.0) {
                flow.position = 0.0; // Clamp to edge start
                
                // Find the next edge connected to the "from" vertex
                const connectedEdges = this.foamEdges.filter((e, idx) => 
                    (e.from === edge.from || e.to === edge.from) && 
                    idx !== flow.edge
                );
                
                if (connectedEdges.length > 0) {
                    // Pick a connected edge to move to
                    const nextEdgeIdx = Math.floor(Math.random() * connectedEdges.length);
                    const nextEdge = connectedEdges[nextEdgeIdx];
                    flow.edge = this.foamEdges.indexOf(nextEdge);
                    
                    // Start at the correct end of the new edge
                    flow.position = nextEdge.from === edge.from ? 0.0 : 1.0;
                    
                    // Make sure the velocity is in the correct direction to move away from junction
                    if ((flow.position === 0.0 && flow.velocity < 0) || 
                        (flow.position === 1.0 && flow.velocity > 0)) {
                        flow.velocity = -flow.velocity;
                    }
                } else {
                    // No connected edges, just reverse direction
                    flow.velocity = -flow.velocity;
                    flow.position = 0.01; // Move slightly forward to avoid getting stuck
                }
            }
        }
    }
    
    /**
     * Helper function to find shared points between two triangles
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
     * Add flow to a specific Delaunay edge
     */
    addFlowToEdge(p1, p2, amount) {
        // Keep track of flow on edges between points
        const edgeKey = `${Math.min(p1, p2)}-${Math.max(p1, p2)}`;
        
        if (!this.edgeFlows) {
            this.edgeFlows = {};
        }
        
        if (!this.edgeFlows[edgeKey]) {
            this.edgeFlows[edgeKey] = 0;
        }
        
        this.edgeFlows[edgeKey] += amount;
    }
    
    /**
     * Update Delaunay points based on edge flows
     */
    updateDelaunayPoints(deltaTime) {
        if (!this.edgeFlows) return false;
        
        let maxMovement = 0;
        
        // Calculate point movements
        const movements = Array(this.points.length / 2).fill().map(() => [0, 0]);
        
        // Apply flows to move Delaunay points
        for (const edgeKey in this.edgeFlows) {
            const [p1, p2] = edgeKey.split('-').map(Number);
            const flow = this.edgeFlows[edgeKey];
            
            if (Math.abs(flow) < 0.01) {
                this.edgeFlows[edgeKey] *= 0.99; // Decay small flows
                continue;
            }
            
            // Get point coordinates
            const x1 = this.points[p1 * 2];
            const y1 = this.points[p1 * 2 + 1];
            const x2 = this.points[p2 * 2];
            const y2 = this.points[p2 * 2 + 1];
            
            // Calculate current distance
            const currentDistance = Math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2);
            
            // Target distance based on flow
            let targetDistance = this.equilibriumDistance; // Base distance
            
            if (flow > this.expansionThreshold) {
                // Expansion for high flow
                targetDistance *= this.expansionFactor;
            } else {
                // Contraction for normal flow
                targetDistance *= this.contractionFactor;
            }
            
            // Move points towards target distance
            const factor = (targetDistance / currentDistance - 1) * 0.1;
            const dx = (x2 - x1) * factor;
            const dy = (y2 - y1) * factor;
            
            // Accumulate movements for both points
            movements[p1][0] -= dx;
            movements[p1][1] -= dy;
            movements[p2][0] += dx;
            movements[p2][1] += dy;
            
            // Track maximum movement
            const maxDelta = Math.max(Math.abs(dx), Math.abs(dy));
            maxMovement = Math.max(maxMovement, maxDelta);
            
            // Decay flow
            this.edgeFlows[edgeKey] *= 0.95;
        }
        
        // Apply movements to points
        for (let i = 0; i < movements.length; i++) {
            this.points[i * 2] += movements[i][0];
            this.points[i * 2 + 1] += movements[i][1];
            
            // Keep points within bounds
            const maxWidth = this.width / 2 - 20;
            const maxHeight = this.height / 2 - 20;
            this.points[i * 2] = Math.max(-maxWidth, Math.min(maxWidth, this.points[i * 2]));
            this.points[i * 2 + 1] = Math.max(-maxHeight, Math.min(maxHeight, this.points[i * 2 + 1]));
        }
        
        return maxMovement > 0.05; // Return true if significant movement occurred
    }
    
    /**
     * Updates the simulation by one step
     */
    update(deltaTime) {
        // Update flows along Voronoi edges
        this.updateFlows(deltaTime);
        
        // Update Delaunay points based on edge flows
        const significantChanges = this.updateDelaunayPoints(deltaTime);
        
        // If points moved significantly, recalculate triangulation and foam
        if (significantChanges) {
            this.triangulate();
            this.createFoam();
        }
    }
    
    /**
     * Return the geometry data for Three.js visualization
     */
    getGeometryData() {
        return {
            points: this.points,
            triangles: Array.from(this.delaunay.triangles),
            centers: this.triangleCenters,
            edges: this.foamEdges,
            flows: this.foamFlows
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