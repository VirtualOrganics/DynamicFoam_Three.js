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
                position: 0, // Start at the 'from' vertex
                velocity: this.flowSpeed, // Always positive to move from 'from' to 'to'
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
     * Updates the simulation by one step
     */
    update(deltaTime) {
        // Update flows along edges
        this.updateFlows(deltaTime);
        
        // Update foam structure based on flows
        this.updateFoamStructure();
    }
    
    /**
     * Update particle flows along edges
     */
    updateFlows(deltaTime) {
        for (let i = 0; i < this.foamFlows.length; i++) {
            const flow = this.foamFlows[i];
            const edge = this.foamEdges[flow.edge];
            
            // Ensure velocity is always positive and matches the flowSpeed
            flow.velocity = this.flowSpeed;
            
            // Update position
            flow.position += flow.velocity * deltaTime / edge.length;
            
            // If reached the 'to' vertex (end of edge), move to the next edge
            if (flow.position >= 1) {
                // Add flow to current edge
                edge.flow += this.flowStrength;
                
                // The vertex we've hit is the 'to' vertex of the current edge
                const vertex = edge.to;
                
                // Find all edges connected to this vertex
                const connectedEdges = this.foamEdges.filter((e, idx) => 
                    (e.from === vertex || e.to === vertex) && idx !== flow.edge
                );
                
                if (connectedEdges.length > 0) {
                    // Pick a random connected edge
                    const nextEdge = connectedEdges[Math.floor(Math.random() * connectedEdges.length)];
                    const nextEdgeIndex = this.foamEdges.indexOf(nextEdge);
                    flow.edge = nextEdgeIndex;
                    
                    // If the new edge has the vertex as its 'to', reverse the direction by starting at 1
                    // Otherwise, start at 0 (forward direction)
                    if (nextEdge.to === vertex) {
                        this.foamEdges[nextEdgeIndex] = {
                            from: nextEdge.to,
                            to: nextEdge.from,
                            flow: nextEdge.flow,
                            length: nextEdge.length
                        };
                        flow.position = 0; // Start at the new 'from'
                    } else {
                        flow.position = 0; // Start at the 'from' vertex of the new edge
                    }
                    
                    // Velocity remains positive
                    flow.velocity = this.flowSpeed;
                } else {
                    // If no connected edges, reset to start of the current edge
                    flow.position = 0;
                    flow.velocity = this.flowSpeed;
                }
            }
        }
    }
    
    /**
     * Update foam structure based on flows
     */
    updateFoamStructure() {
        // Adjust points based on edge flows
        for (const edge of this.foamEdges) {
            // Skip edges with little flow
            if (Math.abs(edge.flow) < 0.01) continue;
            
            const fromCenter = this.triangleCenters[edge.from];
            const toCenter = this.triangleCenters[edge.to];
            
            // Calculate distance between centers
            const currentDistance = this.calculateDistance(fromCenter, toCenter);
            
            // Target distance based on flow
            let targetDistance = this.equilibriumDistance;
            
            if (edge.flow > this.expansionThreshold) {
                // Expansion for high flow
                targetDistance *= this.expansionFactor;
            } else {
                // Contraction for normal flow
                targetDistance *= this.contractionFactor;
            }
            
            // Move points towards target distance
            if (Math.abs(currentDistance - targetDistance) > 0.1) {
                const dx = toCenter[0] - fromCenter[0];
                const dy = toCenter[1] - fromCenter[1];
                const factor = (targetDistance / currentDistance - 1) * 0.1;
                
                // Move the points
                fromCenter[0] -= dx * factor;
                fromCenter[1] -= dy * factor;
                toCenter[0] += dx * factor;
                toCenter[1] += dy * factor;
            }
            
            // Decay flow over time
            edge.flow *= 0.99;
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